from collections import deque
import random, os
import wx
import numpy as np
import epics
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

from datetime import datetime
from scipy.stats.stats import pearsonr

#from epics.wx import MotorPanel

from slic.gui.widgets import LabeledMathEntry 
from slic.gui.widgets import make_filled_hbox, make_filled_vbox, EXPANDING
from slic.gui.widgets.plotting import PlotPanel

from bstrd import BS, bsstream
#from bsread import source

# config
frames_num   =50

# channel
ph_num_name = 'SATES30-RIXS-CAM01:EVENT_NUM'
I_int_name  = 'SATES30-RIXS-CAM01:EVENT_I_INTERP'
J_int_name  = 'SATES30-RIXS-CAM01:EVENT_J_INTERP'


#create channels
ph_num_ch = BS(ph_num_name)
I_int_ch  = BS(I_int_name)
J_int_ch  = BS(J_int_name)


iso_format = "%H:%M:%S"

class MainPanel(wx.Panel):

    def __init__(self, parent):
        super().__init__(parent)

        self.slope = 2.33581e-03
        self.c2 = 2.381609e-6
        self.pb=1
        self.ROI = [[0,1800],[0,1800]]



        #self.evts  = np.empty((nshots,256))
        self.ph = np.empty(frames_num)
        self.phs= deque(maxlen=100)

        self.col_i_tot = deque([],frames_num)
        self.row_i_tot = deque([],frames_num)
        self.row_curv_corr_i_tot = deque([],frames_num)
        self.spc_1d = deque([],1)
        self.x_pxl  = deque([],1)



        self.plot_int     = plot_int = PlotPanel(self, figsize=(3,1))
        self.plot_spectra = plot_spectra = PlotPanel(self, figsize=(3,1))
        self.plot_image   = plot_image = PlotPanel(self, figsize=(3,1))

        plots1 = (plot_int, )
        plots2 = (plot_spectra, )
        plots3 = (plot_image, )

        hb_plot1 = make_filled_hbox(plots1)
        hb_plot2 = make_filled_hbox(plots2)
        hb_plot3 = make_filled_hbox(plots3)

        btn_clearQ = wx.Button(self, label="Clear plots")
        btns = (btn_clearQ,)
        hb_btns = make_filled_hbox(btns)


       
        widgets = (hb_plot1,hb_plot2, hb_plot3, )
        box = make_filled_vbox(widgets, border=1)
        self.SetSizerAndFit(box)

        #btn_clearQ.Bind(wx.EVT_BUTTON, self.on_click_clearQ)

        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_update, self.timer)
        self.timer.Start(100)



    def on_update(self, _event):

        self.time = datetime.now().strftime(iso_format) 
        
        self.ph  = np.empty(frames_num)
        

        for i in range(frames_num):
            event_num=ph_num_ch.get()
            self.ph[i] = event_num
            
            col_i=J_int_ch.get()
            self.col_i_tot.append(col_i[0:event_num])

            row_i=I_int_ch.get()
            self.row_i_tot.append(row_i[0:event_num])

            row_curv_corr = row_i[0:event_num] -self.slope*col_i[0:event_num]  -self.c2*col_i[0:event_num]**2
            self.row_curv_corr_i_tot.append(row_curv_corr)
            

            next(bsstream)
            wx.GetApp().Yield()
        
        self.phs.append(np.mean(self.ph, axis=0))

        pxY_num   = self.ROI[1][1]
        startY_px = self.ROI[1][0]
        endY_px   = self.ROI[1][0]+self.ROI[1][1]-1
 
        row_flat  = np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel()
        spc_1d, b = np.histogram(row_flat, bins=int(pxY_num/self.pb), range=(startY_px, endY_px)) #pb is pixel binning
        
        x_pxl = .5*(b[:-1] + b[1:]) #+ B*px_col/2
        self.spc_1d.append(spc_1d)
        self.x_pxl.append(x_pxl)

        self.draw_plot()
  

    def draw_plot(self):


        self.plot_int.clear()
        self.plot_int.set_xlabel('shots')
        self.plot_int.set_ylabel('photon num')

        self.plot_int.plot(np.arange(0, len(self.phs)), self.phs, '-', color='purple')
        
        self.plot_spectra.clear()
        self.plot_spectra.plot(np.asarray(self.x_pxl)[0], np.asarray(self.spc_1d)[0], '-o'  )
        
        self.plot_image.clear()
        self.plot_image.plot( np.concatenate(np.asarray(self.row_curv_corr_i_tot)).ravel(), np.concatenate(np.asarray(self.col_i_tot)).ravel(),  'o', markersize=.5 ) 


        self.plot_int.draw() 
        self.plot_spectra.draw()
        self.plot_image.draw()  




 


