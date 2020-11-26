import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def minmax(_minmax, nd_):
  _min, _max = min(_minmax), max(_minmax)
  d_ = (_max-_min)*0.01
  _min, _max = _min-d_*nd_, _max+d_*nd_
  d_ = (_max-_min)*0.01
  return _min, _max, d_

def draw_ticks(ax):  
  x_majors = ax.xaxis.get_majorticklocs()
  x_minor  = ((max(x_majors) - min(x_majors)) / float(len(x_majors)-1))/5.0  
  ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=x_minor)) 
  
  y_majors = ax.yaxis.get_majorticklocs()
  y_minor  =  ((max(y_majors) - min(y_majors)) / float(len(y_majors)-1))/5.0  
  ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=y_minor))
  
  ax.tick_params('both', length=4, width=0.8, which='major', direction='in')
  ax.tick_params('both', length=2, width=0.8, which='minor', direction='in')
  ax.xaxis.set_ticks_position('both')
  ax.yaxis.set_ticks_position('both')

def lambda_theta(err_x, err_y, cov_xy):

  cov_xx, cov_yy = err_x**2, err_y**2
  rho = cov_xy/np.sqrt(cov_xx*cov_yy) 
  cov = np.zeros(shape=(2,2))
  cov[0][0] = cov_xx
  cov[0][1] = cov_xy
  cov[1][0] = cov_xy
  cov[1][1] = cov_yy
  
  lambda_, v = np.linalg.eig(cov)
  lambda_ = np.sqrt(lambda_)
  
  theta = np.degrees(np.arctan2(*v[:,0][::-1]))
                                
  return lambda_, theta
  
def plot_ellipse(ax, x, ex, y, ey, cov_xy, color, lw):
  
  n_sigma = 1.0
  lambda_, theta = lambda_theta(ex, ey, cov_xy)
  diameter_hor = 2.0*lambda_[0]
  diameter_ver = 2.0*lambda_[1]
  
  ell = Ellipse(xy=(x, y), width=diameter_hor*n_sigma, height=diameter_ver*n_sigma, angle=theta, edgecolor=color, lw=lw)
  ell.set_facecolor('none')
  ax.add_patch(ell)
  ax.add_artist(ell)

def plot_C3(panels, self, figure_name='', sn='', time_label='Time'):
  
  t_BVI = self.t_BVI
  tmin, tmax = self.t_min, self.t_max
  EhBV_i, err_EhBV_i = self.EhBV_i, self.err_EhBV_i
  ssd_EhBV, EhBV, err_EhBV = self.ssd_EhBV, self.EhBV, self.err_EhBV['stat+sys']
  
  ms=6
  ecolor, elw = 'gray', 0.2
  
  colors  = {'B':'royalblue', 'V':'lime', 'I':'r'}
  symbols = {'B':'o', 'V':'s', 'I':'d'}
  panels = panels.split('+')
  orientation = 'v'
  if 'c3' in panels or 'c3-zoom' in panels:  orientation = 'h'
  
  if len(panels) == 1:
      fig = plt.figure(figsize=(4, 3.5))
      fig.subplots_adjust(left=0.175, bottom=0.11, right=0.95, top=0.985)
  if len(panels) == 2:
      if orientation == 'h':
          ms=6
          scale = 0.5
          fig = plt.figure(figsize=(16.0*scale, 7.0*scale))
          fig.subplots_adjust(wspace=0.26, left=0.088, bottom=0.11, right=0.975, top=0.985)
      if orientation == 'v':
          ms=6
          fig = plt.figure(figsize=(5.0, 5.0))
          fig.subplots_adjust(wspace=0.17, left=0.15, bottom=0.075, right=0.975, top=0.985)
  
  ax = {}
  mask_in, mask_out = (t_BVI['t']>=tmin) & (t_BVI['t']<=tmax), (t_BVI['t']<tmin) | (t_BVI['t']>tmax)
  discarded_points = True
  if True not in mask_out: discarded_points = False
  lc_empty_symbols = False
  for i_panel in range(len(panels)):
      
      panel = panels[i_panel]
      
      if orientation == 'h':  ax[panel] = fig.add_subplot(1,len(panels),i_panel+1)
      if orientation == 'v':  ax[panel] = fig.add_subplot(len(panels),1,i_panel+1)
      
      ts = t_BVI['t']
      
      #light curves
      if panel in ['lc', 'lc-zoom']:
          x_minmax, y_minmax = np.array([]), np.array([])
          l = {}
          for band in 'BVI':
              l[band] = {}
              
              m  = t_BVI[band][mask_in]
              em = t_BVI['e'+band][mask_in]
              t  = ts[mask_in]
                            
              ax[panel].plot([t, t], [m-em, m+em], '-', color=ecolor, lw=elw)
              l[band]['filled'], = ax[panel].plot(t, m, symbols[band], color=colors[band], markeredgecolor='k', ms=ms, mew=0.2, label='$'+band+'$')
              
              xminmax = np.array([min(ts), max(ts)])
              yminmax = np.array([min(t_BVI[band]-t_BVI['e'+band]), max(t_BVI[band]+t_BVI['e'+band])])
              if panel == 'lc-zoom':
                  xminmax = np.array([min(t), max(t)])
                  yminmax = np.array([min(m), max(m)])
                  
              x_minmax = np.append(x_minmax, xminmax)
              y_minmax = np.append(y_minmax, yminmax)
                  
              if discarded_points and panel != 'lc-zoom':
                  lc_empty_symbols = True
                  m  = t_BVI[band][mask_out]
                  em = t_BVI['e'+band][mask_out]
                  t  = ts[mask_out]
                  ax[panel].plot([t, t], [m-em, m+em], '-', color=ecolor, lw=elw)
                  l[band]['empty'], = ax[panel].plot(t, m, symbols[band], color='none', markeredgecolor=colors[band], ms=ms, mew=0.3)
      
      #Eh
      if panel in ['Eh', 'Eh-zoom']:
          
          Eh  =     EhBV_i[mask_in]
          eEh = err_EhBV_i[mask_in]
          t   =         ts[mask_in]
          
          ax[panel].plot([min(t),max(t)], [EhBV, EhBV], '-k')
          ax[panel].plot([min(t),max(t)], [EhBV-ssd_EhBV, EhBV-ssd_EhBV], '--k')
          ax[panel].plot([min(t),max(t)], [EhBV+ssd_EhBV, EhBV+ssd_EhBV], '--k')
          
          ax[panel].plot([t,t], [Eh-eEh, Eh+eEh], '-', color=ecolor, lw=elw)
          ax[panel].plot(t, Eh, 'o', color='lime', markeredgecolor='k', ms=ms, mew=0.2)
          
          x_minmax = ts.copy()
          y_minmax = [min(EhBV_i-err_EhBV_i),max(EhBV_i+err_EhBV_i)]
          if panel == 'Eh-zoom':
              x_minmax = t.copy()
              y_minmax = [min(Eh-eEh),max(Eh+eEh)]
          
          if discarded_points and panel != 'Eh-zoom':
              Eh  =     EhBV_i[mask_out]
              eEh = err_EhBV_i[mask_out]
              t   =         ts[mask_out]
              ax[panel].plot([t,t], [Eh-eEh, Eh+eEh], '-', color=ecolor, lw=elw)
              ax[panel].plot(t, Eh, 'o', color='none', markeredgecolor='gray', ms=ms, mew=0.3)
          
          if i_panel == 0:  ax[panel].plot(-99, 99, ',w', label=sn)
          N_str = str(len(EhBV_i))
          N_text = '$N\!=\!$'+N_str
          
          ssd_EhBV_str = str(round(ssd_EhBV,3))
          ssd_EhBV_text = '$\hat{\sigma}\!=\!$'+ssd_EhBV_str
          ax[panel].plot(-99, 99, ',w', label=N_text+'; '+ssd_EhBV_text)
          
          if False:
              EhBV_str= str(EhBV)
              while len(EhBV_str)<5: EhBV_str = EhBV_str + '0'
              eEhBV_str= str(err_EhBV)
              while len(eEhBV_str)<5: eEhBV_str = eEhBV_str + '0'
              
              EhBV_text = '$E^{\mathrm{h}}_{B\!-\!V}\!=\!'+EhBV_str+'\!\pm\!'+eEhBV_str+'$'
              ax[panel].plot(-99, 99, ',w', label=EhBV_text)
          
          
      #color-color plot
      if panel in ['c3', 'c3-zoom']:
          AK_colors = self.AK_colors          
          
          x_minmax = np.append(AK_colors['BV']-AK_colors['eBV'], AK_colors['BV']+AK_colors['eBV'])
          y_minmax = np.append(AK_colors['VI']-AK_colors['eVI'], AK_colors['VI']+AK_colors['eVI'])
          
          BVs, eBVs, VIs, eVIs, cBV_VIs = AK_colors['BV'][mask_in], AK_colors['eBV'][mask_in], AK_colors['VI'][mask_in], AK_colors['eVI'][mask_in], AK_colors['cBV_VI'][mask_in]
          ax[panel].plot(BVs, VIs, 'o', color='lime', mec='k', ms=ms, mew=0.2)
          for x, ex, y, ey, cov_xy in zip(BVs, eBVs, VIs, eVIs, cBV_VIs):
              plot_ellipse(ax[panel], x, ex, y, ey, cov_xy, 'gray', elw)
          
          x_sl, y_sl = BVs.copy(), VIs.copy()
          
          if discarded_points and panel != 'c3-zoom':
              BVs, eBVs, VIs, eVIs, cBV_VIs = AK_colors['BV'][mask_out], AK_colors['eBV'][mask_out], AK_colors['VI'][mask_out], AK_colors['eVI'][mask_out], AK_colors['cBV_VI'][mask_out]
              ax[panel].plot(BVs, VIs, 'o', color='none', mec='dimgray', ms=ms)
              for x, ex, y, ey, cov_xy in zip(BVs, eBVs, VIs, eVIs, cBV_VIs):
                  plot_ellipse(ax[panel], x, ex, y, ey, cov_xy, 'gray', elw)
                       
          #jackknife
          slope_text = ''
          if len(x_sl) >= 3:
              pars = np.polyfit(x_sl, y_sl, 1)[::-1]
              slope_text = str(round(pars[1],3))
              if len(slope_text)==4: slope_text=slope_text+'0'
              
              slopes = np.array([])
              for i in range(0, len(x_sl)):
                  
                  if i == 0:
                      x_jk, y_jk = x_sl[1:len(x_sl)], y_sl[1:len(y_sl)]
                  elif i == len(x_sl) - 1:
                      x_jk, y_jk = x_sl[0:i], y_sl[0:i]
                  else:
                      x_jk = np.append(x_sl[0:i], x_sl[i+1:len(x_sl)])
                      y_jk = np.append(y_sl[0:i], y_sl[i+1:len(y_sl)])
                  slope = np.polyfit(x_jk, y_jk, 1)[::-1][1]
                  slopes = np.append(slopes, slope)
              
              residuals = slopes - pars[1]
              err_slope = np.sqrt(np.sum(residuals**2)/float(len(slopes)-1))
              err_slope_text = str(round(err_slope,3))
              if len(err_slope_text)==4: err_slope_text=err_slope_text+'0'
              slope_text = slope_text + '$\pm$'+err_slope_text
              
              x_fit = np.array([min(x_sl), max(x_sl)])
              y_fit = pars[0] + pars[1]*x_fit
              ax[panel].plot(x_fit, y_fit, '--k', zorder=0)
          
          
          if i_panel == 0:  ax[panel].plot(-99, 99, ',w', label=sn)
          if slope_text != '':  ax[panel].plot(-99, 99, ',w', label='Slope='+slope_text)
          
          x_fit = np.array([min(x_sl), max(x_sl)])
          y_fit = 0.108 + 0.45*x_fit
          ax[panel].plot(x_fit, y_fit, ':b', zorder=0, label='Reddening-free C3')
          
          if panel == 'c3-zoom':
              x_minmax = x_sl
              y_minmax = np.append(y_sl, y_fit)
      
      xmin, xmax, dx = minmax(x_minmax, 5.0)
      ymin, ymax, dy = minmax(y_minmax, 5.0)
      ax[panel].set_xlim(xmin, xmax)
      ax[panel].set_ylim(ymin, ymax)
      draw_ticks(ax[panel])
      if panel in ['lc', 'lc-zoom']:
          ax[panel].invert_yaxis()  
          ax[panel].set_xlabel(time_label)
          ax[panel].set_ylabel('Apparent magitude')
          mpl.rcParams['legend.handlelength'] = 2.0
          
          if lc_empty_symbols == True :  handles = [(l['B']['filled'],l['B']['empty']), (l['V']['filled'],l['V']['empty']), (l['I']['filled'],l['I']['empty'])]
          if lc_empty_symbols == False:  handles = [l['B']['filled'], l['V']['filled'],l['I']['filled']]
          _, labels = ax[panel].get_legend_handles_labels()
          
          if i_panel == 0:  ax[panel].legend(handles = handles, labels=labels, loc='best', ncol=1, handletextpad=0.3, columnspacing=0.1, edgecolor='none', title=sn, handler_map = {tuple: mpl.legend_handler.HandlerTuple(None)})
          if i_panel == 1:  ax[panel].legend(loc='best', ncol=3, handletextpad=-0.3, columnspacing=0.1, edgecolor='none')
      if panel in ['Eh', 'Eh-zoom']:
          ax[panel].set_xlabel(time_label)
          ax[panel].set_ylabel('$E^{\mathrm{h}}_{B\!-\!V,i}$')
          mpl.rcParams['legend.handlelength'] = 0.5
          ax[panel].legend(loc='best', ncol=1, handletextpad=-0.3, labelspacing=0.06, edgecolor='none')
      if panel in ['c3', 'c3-zoom']:
          ax[panel].set_xlabel('$B\!-\!V$')
          ax[panel].set_ylabel('$V\!-\!I$')
          mpl.rcParams['legend.handlelength'] = 1.5
          ax[panel].legend(loc='best', ncol=1, labelspacing=0.06, borderpad=0.2, handletextpad=0.0, edgecolor='none')
      
  if figure_name != '':  fig.savefig(figure_name, dpi=100)

