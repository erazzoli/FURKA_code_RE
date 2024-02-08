from bstrd import BS, bsstream
from bsread import source
from collections import deque
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from threading import Thread
from time import sleep
from scipy.optimize import curve_fit

#import warnings
#warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 



#Define the Gaussian function
def gauss(x, H, A, x0, sigma):
    return H + A**2 * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
 


#Define the Gaussian function
def gauss_test(x, A ):
    return  A * np.exp(-(x -10) ** 2 / (2 * 1 ** 2))


plt.ion()

 





class SPC_events_num:
    def __init__(self,frames_num = 10):
          self.event_num_tot = deque([],frames_num)
          self.col_i_tot = deque([],1)
          self.row_i_tot = deque([],1)
          self.accumulator = Thread(target=self.run_continuously)
          self.accumulator.start()

    def run_continuously(self):
           bs_pids      = BS("pid") 
           bs_col_i     = BS("SATES30-RIXS-CAM01:EVENT_I_INTERP")
           bs_row_i     = BS("SATES30-RIXS-CAM01:EVENT_J_INTERP")
           bs_event_num = BS("SATES30-RIXS-CAM01:EVENT_NUM") 
           bs_I0        = BS("SATES30-LSCP10-FNS:CH0:VAL_GET")
           bs_FY        = BS("SATES30-LSCP10-FNS:CH2:VAL_GET")
          

           for _ in bsstream:
                ix = bs_pids.value
              
                col_i = bs_col_i.value
                if col_i is None:
                    continue
                
                row_i = bs_row_i.value
                if row_i is None:
                    continue

                event_num = bs_event_num.value
                if event_num is None:
                    continue
               
                I0 = bs_I0.value

                FY = bs_FY.value
                if FY is None:
                    continue
                
                
                
                self.event_num_tot.append(event_num)
                self.col_i_tot.append(col_i)
                self.row_i_tot.append(row_i)
		#self.spec_ppd.append(prof)
		#self.pos.append(np.argmin(np.gradient(self.sig[-1])))
		#print(self.pos[-1])
		# print(f'pumped_id is {m.data.pulse_id}')
		
    def setup_plot(self):
        self.ax0_plot = self.axs[0].plot(np.asarray(self.event_num_tot), '-o')[0]
        #self.ax1 = self.axs[1].plot(np.asarray(self.row_i_tot)[0, 0:10],  np.asarray(self.col_i_tot)[0, 0:10], '-o')[0]
        

    def update_plot(self,dum):
        #self.axs[0].relim()
        #self.axs[0].autoscale_view(True,True,True)
        self.ax0_plot.set_ydata(np.asarray(self.event_num_tot))
        #self.ax1.set_data(np.asarray(self.row_i_tot)[0, 0:10],  np.asarray(self.col_i_tot)[0, 0:10])
        return self.ax0_plot

    def plot_animation(self,name='SPC online Image',animate=True):
        #if len(self.sig)<1:
        #    print('no signals yet')
        #    return
        self.fig,self.axs = plt.subplots(2,1,sharex=True,num=name)
        # self.fig.clf()
        #self.ax = self.fig.add_subplot(111)
        if animate:
            self.ani = FuncAnimation(self.fig,self.update_plot,init_func=self.setup_plot)
            plt.show()





class SPC_Image:
    def __init__(self,frames_num = 5, slope = -0.005169, pb = 1, ROI=[[0,1800],[0,1800]], gauss_fit = "True", fit_start_pix= 680, fit_end_pix=760,  pix_to_Energy = 1000/20 ):
          # pix_to_Energy is in meV  
          self.event_num_tot = deque([],1)
          self.frames_num = frames_num
          self.slope = slope
          self.pb=pb
          self.ROI = ROI
          self.gauss_fit = gauss_fit
          self.fit_start_pix=fit_start_pix
          self.fit_end_pix=fit_end_pix
          self.pix_to_Energy=pix_to_Energy


          self.col_i_tot = deque([],frames_num)
          self.row_i_tot = deque([],frames_num)
          self.row_curv_corr_i_tot = deque([],frames_num)
          self.spc_1d = deque([],1)
          self.x_pxl  = deque([],1)
          self.fit_y  = deque([],1)
          self.fit_x  = deque([],1)
          self.pfit   = deque([],1)
          self.Eres   = deque([],1)
          self.accumulator = Thread(target=self.run_continuously)
          self.accumulator.start()

    def run_continuously(self):
           bs_pids      = BS("pid") 
           bs_col_i     = BS("SATES30-RIXS-CAM01:EVENT_I_INTERP")
           bs_row_i     = BS("SATES30-RIXS-CAM01:EVENT_J_INTERP")
           bs_event_num = BS("SATES30-RIXS-CAM01:EVENT_NUM") 
           bs_I0        = BS("SATES30-LSCP10-FNS:CH0:VAL_GET")
           bs_FY        = BS("SATES30-LSCP10-FNS:CH2:VAL_GET")
          

           for _ in bsstream:
                ix = bs_pids.value
              
                col_i = bs_col_i.value
                if col_i is None:
                    continue
                
                row_i = bs_row_i.value
                if row_i is None:
                    continue

                event_num = bs_event_num.value
                if event_num is None:
                    continue
               
                I0 = bs_I0.value

                FY = bs_FY.value
                if FY is None:
                    continue
                
                
                self.event_num_tot.append(event_num)
                self.col_i_tot.append(col_i[0:event_num])
                self.row_i_tot.append(row_i[0:event_num])
                row_curv_corr = row_i[0:event_num] -self.slope*col_i[0:event_num]
                self.row_curv_corr_i_tot.append(row_curv_corr)
                
                pxY_num   = self.ROI[1][1]
                startY_px = self.ROI[1][0]
                endY_px   = self.ROI[1][0]+self.ROI[1][1]-1

                spc_1d, b = np.histogram(np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel(), bins=int(pxY_num/self.pb), range=(startY_px, endY_px)) #pb is pixel binning
                x_pxl = .5*(b[:-1] + b[1:]) #+ B*px_col/2
                self.spc_1d.append(spc_1d)
                self.x_pxl.append(x_pxl)
                
                if self.gauss_fit=="True":                
                   
                   fit_start_i = int(self.fit_start_pix/self.pb)
                   fit_end_i   = int(self.fit_end_pix/self.pb)     
                   
                   p0=[0,100,x_pxl[spc_1d.argmax()], 5]
                                    
      
                   
                   pfit, covariance = curve_fit(gauss, x_pxl,  spc_1d, p0=p0)

                   self.pfit.append(pfit)
                   self.Eres.append(round( self.pfit[0][3]*2.35*self.pix_to_Energy, 3))           
                   self.fit_y.append(gauss(x_pxl, *pfit))
                   self.fit_x.append(x_pxl)
		

		
    def setup_plot(self):
        #self.ax0_plot = self.axs.plot(np.asarray(self.event_num_tot), '-o')[0]
        #self.ax0_plot = self.axs.plot(np.asarray(self.row_i_tot)[0, 0:10],  np.asarray(self.col_i_tot)[0, 0:10], '-o')[0]
        #self.ax0_plot     = self.axs[0].plot(np.asarray(self.col_i_tot, dtype=object)[0], np.asarray(          self.row_i_tot, dtype=object)[0],  'o', markersize=1 )[0]
        self.ax1_plot     = self.axs[1].plot( np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel(), np.concatenate(np.asarray(self.col_i_tot)).ravel(),  'o', markersize=.5 )[0]
        self.ax2_plot     = self.axs[0].plot( np.asarray(self.x_pxl)[0], np.asarray(self.spc_1d)[0],  '-o', markersize=.5)[0]
        if self.gauss_fit=="True":    
          self.ax2_plot_fit = self.axs[0].plot( np.asarray(self.fit_x)[0],np.asarray(self.fit_y)[0],  '-o', markersize=1 )[0]
          self.ax2_text = self.axs[1].text(250, 1950, "E res: " + str(self.Eres[0]) + " meV")
        self.ax2_text2 = self.axs[1].text(1200,250, " Image is sum of: \n" + str(self.frames_num)+ "frames \n ph per frame =" +  str(self.event_num_tot[0]) + " ph")


        #self.axs[0].set_xlim([0, 500])
        #self.axs[0].set_ylim([0, 200])
        #self.axs[1].set_xlim([0, 500])
        #self.axs[1].set_ylim([0, 200])
        

    def update_plot(self,dum):
 
        #self.axs[0].relim()
        self.axs[0].autoscale_view(True,False,True)     
        #self.ax0_plot.set_ydata(np.asarray(self.event_num_tot))
        #self.ax0_plot.set_data(np.asarray(self.col_i_tot,dtype=object)[0], np.asarray(          self.row_i_tot,dtype=object)[0] )
        self.ax1_plot.set_data( np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel(), np.concatenate(np.asarray(self.col_i_tot)).ravel())
        self.ax2_plot.set_data( np.asarray(self.x_pxl, dtype=object)[0], np.asarray(self.spc_1d, dtype=object)[0])
        if self.gauss_fit=="True":    
          self.ax2_plot_fit.set_data( np.asarray(self.fit_x)[0], np.asarray( self.fit_y)[0])
          self.ax2_text.set_text( "E res: " + str(self.Eres[0]) + " meV")
        self.ax2_text2.set_text( "Image is sum of " + str(self.frames_num)+ "frames \n ph per frame =" +  str(self.event_num_tot[0]) + " ph")


        return self.ax1_plot

    def plot_animation(self,name='SPC online analysis',animate=True):
        #if len(self.sig)<1:
        #    print('no signals yet')
        #    return
        self.fig,self.axs = plt.subplots(2,1,sharex=True,num=name)
        # self.fig.clf()
        # self.ax = self.fig.add_subplot(111)
        if animate:
            self.ani = FuncAnimation(self.fig,self.update_plot,init_func=self.setup_plot)
            plt.show()


SPC_events    = SPC_events_num(frames_num = 20)
#SPC_run_V2     = SPC_Image(frames_num = 10, slope = 0.00651614803557518, pb =2, gauss_fit = "False", fit_start_pix= 750, fit_end_pix=870)

sleep(5)

SPC_events.plot_animation()
#SPC_run_V2.plot_animation()


