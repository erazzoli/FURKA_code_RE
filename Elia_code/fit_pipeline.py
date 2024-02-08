import sys
import numpy as np
from logging import getLogger
_logger = getLogger(__name__)

def process(data, pulse_id, timestamp, params):
     try:
          my_parameter = params["my_parameter"]
          event_ch    = data["SATES30-RIXS-CAM01:EVENT_CHARGE"]
          event_eta_x = data["SATES30-RIXS-CAM01:EVENT_ETA_X"]
          event_eta_y = data["SATES30-RIXS-CAM01:EVENT_ETA_Y
          event_I     = data["SATES30-RIXS-CAM01:EVENT_I
          event_I_INT = data["SATES30-RIXS-CAM01:EVENT_I_INTERP
          event_J     = data["SATES30-RIXS-CAM01:EVENT_J
          event_J_INT = data["SATES30-RIXS-CAM01:EVENT_J_INTERP
          event_num   = data["SATES30-RIXS-CAM01:EVENT_NUM


          frames_num= 10
          event_I_tot=np.array([])
          event_J_tot=np.array([])
          event_charge_tot=np.array([])
          event_etax_tot=np.array([])
          event_etay_tot=np.array([])
          event_I_interp_tot=np.array([])
          event_J_interp_tot=np.array([])
    
          event_num_tot=0

          row_curv_corr = np.array([])


    # add together 
    frames_num = subset["SATES30-RIXS-CAM01:EVENT_NUM"].pids.shape[0]
    for i in range(frames_num):
        event_num_i= subset["SATES30-RIXS-CAM01:EVENT_NUM"].data[i]
        event_I_tot= np.append(event_I_tot, subset["SATES30-RIXS-CAM01:EVENT_I_INTERP"].data[i][0:event_num_i])
        event_J_tot= np.append(event_J_tot, subset["SATES30-RIXS-CAM01:EVENT_J_INTERP"].data[i][0:event_num_i])
        event_num_tot+=subset["SATES30-RIXS-CAM01:EVENT_NUM"].data[i]




          data["test"] = 4.0
          data["par"]=my_parameter
           
          

          return data
     except:
         _logger.exception(str(sys.exc_info()[1]))
     
