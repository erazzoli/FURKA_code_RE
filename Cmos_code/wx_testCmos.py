import wx
import numpy as np
import matplotlib.figure as mfigure
import matplotlib.animation as manim
from collections import deque
from threading import Thread
from bsread import source
import numpy as np
from scipy.optimize import curve_fit

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg


#Define the Gaussian function
def gauss(x, H, A, x0, sigma):
    return H + A**2 * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
 

class MyFrame(wx.Frame):
    def __init__(self,frames_num = 10, slope = 2.33581e-03, c2 = 2.381609e-6 , pb = 1, ROI=[[0,1800],[0,1800]], gauss_fit = "False", fit_start_pix= 680, fit_end_pix=760,  pix_to_Energy = 1000/23 ):
        super(MyFrame,self).__init__(None, wx.ID_ANY, size=(800, 600))
        self.fig = mfigure.Figure()
        self.ax = self.fig.add_subplot(111)
        self.canv = FigureCanvasWxAgg(self, wx.ID_ANY, self.fig)
        self.values = []
        self.animator = manim.FuncAnimation(self.fig,self.anim, interval=1000)
        #self.animator = manim.FuncAnimation(self.fig,self.update_plot,init_func=self.setup_plot)
        self.event_num_tot = deque([],frames_num)
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
        self.x_pxl  = deque([],1)
        self.fit_y  = deque([],1)
        self.fit_x  = deque([],1)
        self.pfit   = deque([],1)
        self.Eres   = deque([],1)
        self.accumulator = Thread(target=self.run_continuously)
        self.accumulator.start()







    def run_continuously(self):
        with source(channels=['SATES30-RIXS-CAM01:EVENT_I_INTERP', 'SATES30-RIXS-CAM01:EVENT_J_INTERP','SATES30-RIXS-CAM01:EVENT_NUM'],  dispatcher_url='http://sf-daqsync-15.psi.ch:9100') as s:
            while True:
                m = s.receive()
                ix = m.data.pulse_id
              
                row_i = m.data.data['SATES30-RIXS-CAM01:EVENT_I_INTERP'].value
                if row_i is None:
                    continue
                
                col_i = m.data.data['SATES30-RIXS-CAM01:EVENT_J_INTERP'].value
                if col_i is None:
                    continue

                event_num = m.data.data['SATES30-RIXS-CAM01:EVENT_NUM'].value
                if event_num is None:
                    continue
                
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


 



    def anim(self,i):
        if i%10 == 0:
            self.values = []
        else:
            #self.values.append(self.event_num_tot[-1])
            self.values=self.spc_1d[0]
#        self.ax.clear()
#        self.ax.set_xlim([0,10])
#        self.ax.set_ylim([0,1])        
        #return self.ax.plot(np.arange(1,i%10+1),self.values,'d-')
        data = np.asarray(self.spc_1d)[0]
        return self.ax.plot( data ,'d-')
        #return self.ax.plot( np.asarray(self.x_pxl)[0], np.asarray(self.spc_1d)[0],  '-o', markersize=.5)

wxa = wx.PySimpleApp()
w = MyFrame()
w.Show(True)
wxa.MainLoop()


