from bsread import source
from collections import deque
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from threading import Thread
from time import sleep
import epics
from epics import caput
plt.ion()






def _interpolate_row(y_known, x_known, x_interp):
    y_interp = np.interp(x_interp, x_known, y_known)
    return y_interp



def find_edge(data, step_length=100, edge_type="rising", refinement=1):
    # refine data
    data_length = data.shape[0]
    refined_data = np.apply_along_axis(
        _interpolate_row,
        axis=0,
        arr=data,
        x_known=np.arange(data_length),
        x_interp=np.arange(0, data_length - 1, refinement),
    )

    # prepare a step function and refine it
    step_waveform = np.ones(shape=(step_length,))
    if edge_type == "rising":
        step_waveform[: int(step_length / 2)] = -1
    elif edge_type == "falling":
        step_waveform[int(step_length / 2) :] = -1

    step_waveform = np.interp(
        x=np.arange(0, step_length - 1, refinement), xp=np.arange(step_length), fp=step_waveform
    )

    # find edges
    xcorr = np.apply_along_axis(np.correlate, 0, refined_data, v=step_waveform, mode="valid")
    edge_position = np.argmax(xcorr, axis=0).astype(float) * refinement
    xcorr_amplitude = np.amax(xcorr, axis=0)

    # correct edge_position for step_length
    edge_position += np.floor(step_length / 2)

    return {"edge_pos": edge_position, "xcorr": xcorr, "xcorr_ampl": xcorr_amplitude, "signal":data}








class TtProcessor:
    def __init__(self,Nbg = 10, pix_min = 0,  pix_max = 2000, fb = False,  n_std_filter_I0=0, WL_delay_stage_pos = 26.116,  poly_coef = [-1.42139290e-14,  2.16712147e-10, -8.44448276e-07,  1.57075480e-03, 2.50742803e+01], edge_pos_num = 500  ):
          self.fb = fb
          self.bg = deque([],Nbg)
          self.sig = deque([],1)
          self.pos = deque([],1)
          self.I0 = deque([],1000)
          self.edge_pos = deque([],edge_pos_num)
          self.edge_pos_pids = deque([],edge_pos_num)
          self.edge_pos_I0   = deque([],edge_pos_num)
          self.edge_pos_num=edge_pos_num
          self.dp_corr  = deque([],edge_pos_num) #delay stage equivalent position
          self.spec_ppd = deque([],1)
          self.accumulator = Thread(target=self.run_continuously)
          self.accumulator.start()
          self.pix_min = pix_min
          self.pix_max = pix_max
          self.poly_coef=poly_coef
          self.poly_fit = np.poly1d(self.poly_coef)
          self.n_std_filter_I0 = n_std_filter_I0
          self.STAGEpv = epics.PV("SLAAT31-LMOT-M816:MOT")
          self.counter=0
          self.WL_delay_stage_pos=WL_delay_stage_pos

    def run_continuously(self):
        with source(channels=['SATES31-CAMS187-RIXS1:x_profile','SAT-CVME-TIFALL5:EvtSet', 'SATES30-LSCP10-FNS:CH0:VAL_GET' ]) as s:
            while True:
                
                
                m = s.receive()
                ix = m.data.pulse_id
              
                prof = m.data.data['SATES31-CAMS187-RIXS1:x_profile'].value
                if prof is None:
                    continue

                codes = m.data.data['SAT-CVME-TIFALL5:EvtSet'].value
                if codes is None:
                    continue
                is_reference = codes[154]==1
                try:
                    if (lastgoodix-ix)>1:
                        print(f'missed  {lastgoodix-ix-1} events!')
                except:
                    pass
                lastgoodix = ix
                if is_reference:
                    self.bg.append(prof)
                else:
                    self.sig.append((prof / np.asarray(self.bg).mean(axis=0)))
                    self.spec_ppd.append(prof)
                    self.pos.append(np.argmin(np.gradient(self.sig[-1])))
                    #print(self.pos[-1])
                    # print(f'pumped_id is {m.data.pulse_id}')
                
		
                if(self.fb == True):

                                self.I0.append( m.data.data['SATES30-LSCP10-FNS:CH0:VAL_GET'].value)
                                I0_filter = np.mean(self.I0)+self.n_std_filter_I0*np.std(self.I0)
                                if    (self.I0[-1] >  I0_filter) and (self.I0[-1] >  10000):
                                      self.counter+=1
                                      sig_i=self.sig[-1]
                                      edge_pos_i= find_edge(sig_i[ self.pix_min : self.pix_max ])["edge_pos"] + self.pix_min
                                      self.edge_pos.append( edge_pos_i )

				      
                                      self.dp_corr.append(  self.poly_fit(edge_pos_i))
                                      self.edge_pos_pids.append(ix)
                                      self.edge_pos_I0.append(m.data.data['SATES30-LSCP10-FNS:CH0:VAL_GET'].value)
                                      #print(self.edge_pos[-1])
				
                                dp_corr_avg = np.mean(self.dp_corr)
                                WL_delay_stage_pos = 18.90
                                #print(len(self.dp_corr))
                                #print(dp_corr_avg)
                                actual_pos= self.STAGEpv.get()

                                if  (len(self.dp_corr)==self.edge_pos_num) and  (abs(dp_corr_avg - self.WL_delay_stage_pos) > 0.005) and (abs(dp_corr_avg - self.WL_delay_stage_pos) < .3) and (np.mean(self.I0)>10000) and ((self.counter>1000)):
                                      print( dp_corr_avg )
                                      print( "Move to "  + str( actual_pos - (dp_corr_avg -self.WL_delay_stage_pos)  ))
                                      self.counter=0
                                      self.STAGEpv.put( actual_pos - (dp_corr_avg - self.WL_delay_stage_pos) )
				      
				      
                      
                      
			



 
    def setup_plot(self):
        self.lh_sig = self.axs[1].plot(self.sig[-1])[0]
        #self.de_sig = self.axs[1].plot(1+np.gradient(self.sig[-2]))[0]
        self.lh_bg = self.axs[0].plot(np.asarray(self.bg).mean(axis=0))[0]
        self.lh_bg_last = self.axs[0].plot(self.bg[-1])[0]
        self.lh_sig_last = self.axs[0].plot(self.spec_ppd[-1])[0]
        if(self.fb == True):  
                self.ax1_text = self.axs[0].text( 600, 30,  " edge pos: " + str(self.edge_pos[-1])  + " I0 = " + str(  np.round(self.edge_pos_I0[-1]/1000000, decimals = 2) )+"e6"  + " dp_corr " + str(  np.round(self.dp_corr[-1], decimals = 3) ) )
                self.ax2_text = self.axs[1].text( 50, 0.95,  " dp avg: " + str(np.mean(self.dp_corr)))

    def update_plot(self,dum):
        self.axs[1].relim()
        self.axs[1].autoscale_view(True,False,True)     
        self.lh_sig.set_ydata(self.sig[-1])
        #self.de_sig.set_ydata(1+np.gradient(self.sig[-1]))
        self.lh_bg.set_ydata(np.asarray(self.bg).mean(axis=0))
        self.lh_bg_last.set_ydata(self.bg[-1])
        self.lh_sig_last.set_ydata(self.spec_ppd[-1])
        if(self.fb == True):
                self.ax1_text.set_text( " edge pos: " + str(self.edge_pos[-1])  + " I0 = " + str(  np.round(self.edge_pos_I0[-1]/1000000, decimals = 2) )+"e6"  + " dp_corr" + str( np.round( self.dp_corr[-1], decimals = 3) ) )
                self.ax2_text.set_text( " dp avg: " + str(np.mean(self.dp_corr)))
        return self.lh_sig

    def plot_animation(self,name='TT online analisys',animate=True):
        if len(self.sig)<1:
            print('no signals yet')
            return
        self.fig,self.axs = plt.subplots(2,1,sharex=True,num=name)
        # self.fig.clf()
        # self.ax = self.fig.add_subplot(111)
        if animate:
            self.ani = FuncAnimation(self.fig,self.update_plot,init_func=self.setup_plot)
            plt.show()



            
tt = TtProcessor()
sleep(5)
tt.plot_animation()






#tt=TtProcessor()
#tt.plot_animation()
