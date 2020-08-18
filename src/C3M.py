import numpy as np
from plot_C3M import plot_C3

def EhBV_estimates(t, BV, VI, err_V, err_BV, err_VI, EGBV, z_helio, t_min=0, t_max=2.5e6): 
  
  #data within [t_min, t_max]
  BV     =     BV[t>=t_min]
  VI     =     VI[t>=t_min]
  err_V  =  err_V[t>=t_min]
  err_BV = err_BV[t>=t_min]
  err_VI = err_VI[t>=t_min]
  t      =      t[t>=t_min]
  BV     =     BV[t<=t_max]
  VI     =     VI[t<=t_max]
  err_V  =  err_V[t<=t_max]
  err_BV = err_BV[t<=t_max]
  err_VI = err_VI[t<=t_max]
  t      =      t[t<=t_max]
  
  #B-V corrected by EG(B-V) and K-correction
  A, B, C = 0.032*EGBV, 0.055*EGBV-1.0, BV-0.994*EGBV
  BV_A  = (-B - np.sqrt(B**2-4.0*A*C))/(2.0*A)
  BV_AK = (BV_A-0.643*z_helio)/(1.0+3.027*z_helio)
  
  #V-I corrected by EG(B-V) and K-correction
  VI_A  = VI - EGBV*(1.368-0.038*BV_A-0.004*BV_A**2)
  VI_AK = VI_A - z_helio*(-0.391+1.646*BV_AK)
  
  #mean B-V color
  mean_BV_AK = np.sum(BV_AK/err_BV**2)/np.sum(1.0/err_BV**2)
  
  #Eh(B-V) estimates
  EhBV_i     = 1.081*VI_AK -0.486*BV_AK - 0.116
  err_EhBV_i = np.sqrt(1.168*err_VI**2+0.237*err_BV**2-1.051*err_V**2)
  
  #weighted average
  w_i = 1.0/err_EhBV_i**2
  EhBV = np.sum(EhBV_i*w_i)/np.sum(w_i)
  
  #sample standard deviation
  residuals = EhBV_i-EhBV
  ssd_EhBV  = np.sqrt(np.sum(residuals**2)/float(len(residuals)-1))
  
  return t, EhBV_i, err_EhBV_i, ssd_EhBV, mean_BV_AK

class Color_Color_Curve_Method:
  def __init__(self, t_BVI, EGBV=0.0, z_helio=0.0, t_min=0.0, t_max=2.5e6):
   
    if z_helio > 0.04:  print('Warning: z_helio is higher than 0.04')
   
    err_ZP_C3 = 0.053
    
    t  = t_BVI['t']
    BV = t_BVI['B']-t_BVI['V']
    VI = t_BVI['V']-t_BVI['I']
    err_V  = t_BVI['eV']
    err_BV = np.sqrt(t_BVI['eB']**2+t_BVI['eV']**2)
    err_VI = np.sqrt(t_BVI['eV']**2+t_BVI['eI']**2)
    
    #quantities to build the Eh plot
    t_EhBV_i_all, EhBV_i_all, err_EhBV_i_all = EhBV_estimates(t, BV, VI, err_V, err_BV, err_VI, EGBV, z_helio)[0:3]
    
    #quantities to compute Eh(B-V)
    EhBV_i, err_EhBV_i, ssd_EhBV, mean_BV_AK = EhBV_estimates(t, BV, VI, err_V, err_BV, err_VI, EGBV, z_helio, t_min=t_min, t_max=t_max)[1:]
    
    #weighted average
    w = 1.0/err_EhBV_i**2
    EhBV = np.sum(EhBV_i*w)/np.sum(w) 
    
    #errors
    err_EhBV_stat = np.sqrt(ssd_EhBV**2+0.006*(mean_BV_AK-0.927*EhBV)**2+0.758*z_helio**2)
    err_EhBV = {}
    err_EhBV['stat']     = round(err_EhBV_stat,3)
    err_EhBV['stat+sys'] = round(np.sqrt(err_EhBV_stat**2+err_ZP_C3**2),3)
    
    self.t_BVI      = t_BVI
    self.t_EhBV_i   = t_EhBV_i_all
    self.EhBV_i     = EhBV_i_all
    self.err_EhBV_i = err_EhBV_i_all
    self.t_min      = t_min
    self.t_max      = t_max
    self.EhBV       = round(EhBV,3)
    self.err_EhBV   = err_EhBV
    self.ssd_EhBV   = ssd_EhBV
      
  def plot(self, panels, sn='', time_label='Time', figure_name=''):
    plot_C3(panels, self, figure_name=figure_name, sn=sn, time_label=time_label)
