from bsread import source
from collections import deque
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from threading import Thread
from time import sleep
from scipy.optimize import curve_fit

import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 



#Define the Gaussian function
def gauss(x, H, A, x0, sigma):
    return H + A**2 * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
 


#Define the Gaussian function
def gauss_test(x, A ):
    return  A * np.exp(-(x -10) ** 2 / (2 * 1 ** 2))


plt.ion()

 



class SPC_Image:
    def __init__(self,frames_num = 200, slope = 2.33581e-03, c2 = 2.381609e-6 , pb = 1, ROI=[[0,1800],[0,300]],
gauss_fit = "True", fit_start_pix= 680, fit_end_pix=760,  pix_to_Energy = 1000/23 ):
          # pix_to_Energy is in meV  
          self.event_num_tot = deque([],frames_num)
          self.pids_SASE = deque([], frames_num)
          self.pids_RIXS = deque([], frames_num)
          self.frames_num = frames_num
          self.slope = slope
          self.c2 = c2
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
          self.spc_1d_tot= deque([],frames_num)
          self.sase_spec_tot= deque([],frames_num)
          self.sase_x = deque([],1)
          self.cc_line_x = deque([],1)
          self.cc_line_y = deque([],1)
          self.x_pxl  = deque([],1)
          self.fit_y  = deque([],1)
          self.fit_x  = deque([],1)
          self.cc     = deque([],1)
          self.pfit   = deque([],1)
          self.Eres   = deque([],1)
          self.accumulator = Thread(target=self.run_continuously)
          self.accumulator.start()

    def run_continuously(self):
        with source(channels=['SATES30-RIXS-CAM01:EVENT_I_INTERP', 'SATES30-RIXS-CAM01:EVENT_J_INTERP','SATES30-RIXS-CAM01:EVENT_NUM','SATOP11-PSAS079:SPECTRUM_X', 'SATOP11-PSAS079:SPECTRUM_Y']) as s:
            while True:
                m = s.receive()
                pid = m.data.pulse_id
              
                
                row_i = m.data.data['SATES30-RIXS-CAM01:EVENT_I_INTERP'].value
                if row_i is None:
                    continue
                
                col_i = m.data.data['SATES30-RIXS-CAM01:EVENT_J_INTERP'].value
                if col_i is None:
                    continue

                event_num = m.data.data['SATES30-RIXS-CAM01:EVENT_NUM'].value
                if event_num is None:
                    continue
                else:
                    self.pids_RIXS.append(pid)



                sase_spec = m.data.data['SATOP11-PSAS079:SPECTRUM_Y'].value
                if sase_spec is None:
                     continue
                else:
                     self.pids_SASE.append(pid)
                
                sase_x = m.data.data['SATOP11-PSAS079:SPECTRUM_X'].value
                if sase_x is None:
                      continue
            


                #print(sase_spec)
                self.sase_spec_tot.append(sase_spec)
                self.sase_x.append(sase_x)                
                self.event_num_tot.append(event_num)
                self.col_i_tot.append(col_i[0:event_num])
                self.row_i_tot.append(row_i[0:event_num])
                row_curv_corr = row_i[0:event_num] -self.slope*col_i[0:event_num]  -self.c2*col_i[0:event_num]**2
                self.row_curv_corr_i_tot.append(row_curv_corr)
                
                pxY_num   = self.ROI[1][1]
                startY_px = self.ROI[1][0]
                endY_px   = self.ROI[1][0]+self.ROI[1][1]-1
                
                           
                                              

                spc_1d, b = np.histogram(np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel(), bins=int(pxY_num/self.pb), range=(startY_px, endY_px)) #pb is pixel binning
                x_pxl = .5*(b[:-1] + b[1:]) #+ B*px_col/2
                self.spc_1d.append(spc_1d)
                self.spc_1d_tot.append(spc_1d)
                self.x_pxl.append(x_pxl)
               
                #shift SASE and RIXS PIDS
                sase_pids = np.array(self.pids_SASE) 
                sase_pids+=1        
                rixs_pids = np.array(self.pids_SASE)
                rixs_pids+=0
                #print(sase_pids)
                #print(rixs_pids) 

                #manual drop missing with delta
                _pids, ind_SASE, ind_RIXS = np.intersect1d(sase_pids, rixs_pids, return_indices=True)
                print(ind_RIXS)
                #print(ind_SASE)
                #print(ind_RIXS)  
                
                #cross correlation here
                RIXS = np.transpose(self.spc_1d_tot)
                SASE = np.transpose(self.sase_spec_tot)
                
                RIXS_dm = RIXS[ind_RIXS]
                SASE_dm = SASE[ind_SASE]

                #self.cc.append(np.cov(RIXS_dm, SASE_dm)) 
                self.cc.append(np.cov(RIXS,SASE)) 
                #np.cov(np.transpose(self.spc_1d_tot), np.transpose(self.sase_spec_tot))
                
                self.cc_line_x.append(np.linspace(0,1999, 2000)) 
                self.cc_line_y.append(self.cc[0][1500,0:2000])                

                p0=[0,1,self.cc_line_x[0][self.cc_line_y[0].argmax()], 10]
                
                #fit, covariance = curve_fit(gauss, self.cc_line_x[0], self.cc_line_y[0], p0=p0)
                #self.pfit.append(pfit)
                #self.Eres.append(round( self.pfit[0][3]*2.35, 3))
                #self.fit_y.append(gauss(self.cc_line_x[0], *pfit))
                #self.fit_x.append(self.cc_line_x[0]) 
 
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
        self.ax3_plot     = self.axs[2].plot( self.sase_x[0],self.sase_spec_tot[0])[0]
        self.ax4_plot     = self.axs[3].plot( self.cc_line_x[0], self.cc_line_y[0] )[0]
        #self.ax4_plot     = self.axs[0].plot( np.asarray(self.fit_x)[0],np.asarray(self.fit_y)[0],  '-o', markersize=1 )[0]

        #self.ax3_plot     = self.axs[2].plot(np.asarray(self.event_num_tot), '-o')[0]
        if self.gauss_fit=="True":    
          self.ax2_plot_fit = self.axs[0].plot( np.asarray(self.fit_x)[0],np.asarray(self.fit_y)[0],  '-o', markersize=1 )[0]
          self.ax2_text = self.axs[1].text(250, 1950, "E res: " + str(self.Eres[0]) + " meV")
        #self.ax2_text2 = self.axs[1].text(1200,250, " Image is sum of: \n" + str(self.frames_num)+ "frames \n ph per frame =" +  str(self.event_num_tot[0]) + " ph")
        


        #self.axs[0].set_xlim([0, 1800])
        #self.axs[0].set_ylim([0, 200])
        #self.axs[1].set_xlim([0, 1800])
        #self.axs[1].set_ylim([0, 1800])
        #self.axs[1].set_ylim([0, 200])
        

    def update_plot(self,dum):
 
        self.axs[0].relim()
        self.axs[0].autoscale_view(True,False,True)    

        #self.axs[1].relim()
        #self.axs[1].autoscale_view(True,False,True)    

        #self.axs[2].relim()
        #self.axs[2].autoscale_view(True,True,True)     
 

        #self.ax0_plot.set_ydata(np.asarray(self.event_num_tot))
        #self.ax0_plot.set_data(np.asarray(self.col_i_tot,dtype=object)[0], np.asarray(          self.row_i_tot,dtype=object)[0] )
        

        self.ax1_plot.set_data( np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel(), np.concatenate(np.asarray(self.col_i_tot)).ravel())
        self.ax2_plot.set_data( np.asarray(self.x_pxl, dtype=object)[0], np.asarray(self.spc_1d, dtype=object)[0])
        #x = np.linspace(0,189,190)
        self.ax3_plot.set_data(self.sase_x[0], self.sase_spec_tot[0])
        self.ax4_plot.set_data(self.cc_line_x[0], self.cc_line_y[0])
        #self.ax4_plot_fit.set_data( np.asarray(self.fit_x)[0], np.asarray( self.fit_y)[0])
        #self.ax4_text.set_text( "E res: " + str(self.Eres[0]) + " meV, Ec[px]:" + str(round(self.pfit[0][2],1)))


        #self.ax3_plot.set_ydata(np.asarray(self.event_num_tot))
        if self.gauss_fit=="True":    
          self.ax2_plot_fit.set_data( np.asarray(self.fit_x)[0], np.asarray( self.fit_y)[0])
          self.ax2_text.set_text( "E res: " + str(self.Eres[0]) + " meV, Ec[px]:" + str(round(self.pfit[0][2],1)))
        #self.ax2_text2.set_text( "Image is sum of " + str(self.frames_num)+ "frames \n ph per frame =" +  str(self.event_num_tot[0]) + " ph")
        


        return self.ax1_plot

    def plot_animation(self,name='SPC online analysis',animate=True):
        #if len(self.sig)<1:
        #    print('no signals yet')
        #    return
        self.fig,self.axs   = plt.subplots(4,1,sharex=True,num=name)
        #self.fig,self.axs   = plt.subplots(3,1, num=name)

        # self.fig.clf()
        # self.ax = self.fig.add_subplot(111)
        if animate:
            self.ani = FuncAnimation(self.fig,self.update_plot,init_func=self.setup_plot)
            plt.show()


SPC_run = SPC_Image(frames_num =5, gauss_fit="False", fit_start_pix= 750, fit_end_pix=870)

sleep(5)

SPC_run.plot_animation()


