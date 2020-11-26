import numpy as np
from plot_C3M import plot_C3

class Color_Color_Curve_Method:
  def __init__(self, t_BVI, EGBV=0.0, z_helio=0.0, t_min=0.0, t_max=2.5e6):
      
    R_V, R_I = 3.1, 1.72
    m_C3 , err_m_C3  =  0.45 , 0.07  #V-I versus B-V C3 slope
    ZP_C3, err_ZP_C3 = -0.116, 0.024 #reddening zero-point for the C3 method
    mean_BV_ZP       =  1.065        #mean B-V color of the SNe used to compute the reddening zero-point
    
    #apparent B-V and V-I colors
    BV, err_BV = t_BVI['B']-t_BVI['V'], np.sqrt(t_BVI['eB']**2+t_BVI['eV']**2)
    VI, err_VI = t_BVI['V']-t_BVI['I'], np.sqrt(t_BVI['eV']**2+t_BVI['eI']**2)
        
    #A_i estimates
    A_i_all, err_A_i_all, BV_AK_i_all, VI_AK_i_all= A_estimates(BV, err_BV, VI, err_VI, t_BVI['eV'], EGBV, z_helio, m_C3, R_V, R_I)
    EhBV_i_all = A_i_all + ZP_C3

    #data within the range [t_min, t_max]
    mask     = (t_BVI['t']>=t_min) & (t_BVI['t']<=t_max)
    A_i      =     A_i_all[mask]
    err_A_i  = err_A_i_all[mask]
    BV_AK_i  = BV_AK_i_all[mask]
    err_BV_i =      err_BV[mask]
 
    #A value that maximizes the log-likelihood of a constant-only model
    if len(A_i) > 1:
        A, ssd_A = likelihood_maximization(A_i, err_A_i, n_pars=1)
    else:
        print('Warning: unable to compute the sample standard deviation (only one point available)')
        A, ssd_A  = A_i[0], err_A_i[0]
    
    #final Eh(B-V)
    EhBV = A + ZP_C3
    
    #error in Eh(B-V) induced by the m_C3 error
    mean_BV_AK = np.sum(BV_AK_i/err_BV_i**2)/np.sum(1.0/err_BV_i**2)
    err_m = np.abs((EhBV-mean_BV_AK+mean_BV_ZP)/(R_V-R_I-m_C3))*err_m_C3
    
    #statistical and systematic error in Eh(B-V)
    err_EhBV_stat = np.sqrt(ssd_A**2+err_m**2)
    err_EhBV_sys  = err_ZP_C3
    err_EhBV      = np.sqrt(err_EhBV_stat**2+err_EhBV_sys**2)
    
    self.t_BVI      = t_BVI
    self.EhBV_i     = EhBV_i_all
    self.err_EhBV_i = err_A_i_all
    self.t_min      = t_min
    self.t_max      = t_max
    self.EhBV       = round(EhBV,3)
    self.err_EhBV   = {'stat': round(err_EhBV_stat,3), 'sys':round(err_EhBV_sys,3), 'stat+sys':round(err_EhBV,3) }
    self.ssd_EhBV   = ssd_A
    self.AK_colors  = {'BV': BV_AK_i_all, 'eBV':err_BV, 'VI': VI_AK_i_all, 'eVI':err_VI, 'cBV_VI':-t_BVI['eV']**2}
      
  def plot(self, panels, sn='', time_label='Time', figure_name=''):
    plot_C3(panels, self, figure_name=figure_name, sn=sn, time_label=time_label)

def AG_correction(BV, VI, EGBV):
  
  if EGBV == 0:
      BV_A, VI_A = BV.copy(), VI.copy()
  else:
      A, B, C = 0.032*EGBV, 0.055*EGBV-1.0, BV-0.994*EGBV
      BV_A = (-B - np.sqrt(B**2-4.0*A*C))/(2.0*A)
      VI_A = VI - EGBV*(1.368-0.038*BV_A-0.004*BV_A**2)
  
  return BV_A, VI_A

def K_correction(BV_A, VI_A, z_helio):
  
  if z_helio > 0.04:  
      print('Warning: implemented K-correction is valid for z_helio<0.04 (input z_helio='+str(z_helio)+')')
      print('         Consider providing K-corrected BVI photometry to the C3M.')
    
  BV_AK = (BV_A-0.643*z_helio)/(1.0+3.027*z_helio)
  VI_AK = VI_A - z_helio*(-0.391+1.646*BV_AK)
  
  return BV_AK, VI_AK

def likelihood_maximization(A_i, err_A_i, n_pars=1):

  N, A = len(A_i), np.mean(A_i)
  ssd  = np.sqrt(np.sum((A_i-A)**2)/float(N-n_pars))
  if N == 2:
      e0s = [0.0]
  else:
      e0s = np.linspace(0.0, 2.0*ssd, 201)
      
  m2lnL_min = 1.e99
  for e0 in e0s:
      Var = err_A_i**2+e0**2
      A = np.sum(A_i/Var)/np.sum(1.0/Var)
      m2lnL = np.sum(np.log(Var)+(A_i-A)**2/Var)
      if m2lnL < m2lnL_min:
          m2lnL_min = m2lnL
          A_min     = A
          e0_min    = e0
  if e0_min != 0.0:  n_pars = n_pars + 1
  
  A   = A_min
  ssd = np.sqrt(np.sum((A_i-A)**2)/float(N-n_pars))
  
  return A, ssd

def A_estimates(BV, err_BV, VI, err_VI, err_V, EGBV, z_helio, m_C3, R_V, R_I): 
    
  #(B-V) and (V-I) corrected for EG(B-V) and K-correction
  BV_A, VI_A = AG_correction(BV, VI, EGBV)
  
  #(B-V)_A and (V-I)_A corrected for K-correction
  BV_AK_i, VI_AK_i = K_correction(BV_A, VI_A, z_helio)
  
  #A_i estimates
  den     = R_V-R_I-m_C3
  A_i     = (VI_AK_i - m_C3*BV_AK_i)/den
  err_A_i = np.sqrt((err_VI/den)**2+(err_BV*m_C3/den)**2-2.0*m_C3*(err_V/den)**2)
  
  return A_i, err_A_i, BV_AK_i, VI_AK_i
