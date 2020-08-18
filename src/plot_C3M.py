import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

mpl.rc('axes', labelsize=18.0)
mpl.rcParams.update({'font.size': 18})

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
  
  ax.tick_params('both', length=4, width=1, which='major', direction='in')
  ax.tick_params('both', length=2, width=1, which='minor', direction='in')
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
  t_EhBV_i, EhBV_i, err_EhBV_i = self.t_EhBV_i, self.EhBV_i, self.err_EhBV_i
  ssd_EhBV, EhBV, err_EhBV = self.ssd_EhBV, self.EhBV, self.err_EhBV['stat+sys']
  
  ms=8
  ecolor, elw = 'gray', 0.5
  
  colors  = {'B':'royalblue', 'V':'lime', 'I':'r'}
  symbols = {'B':'o', 'V':'s', 'I':'d'}
  panels = panels.split('+')
  
  if len(panels) == 1:
      fig = plt.figure(figsize=(8, 7))
      fig.subplots_adjust(left=0.14, bottom=0.087, right=0.95, top=0.985)
  if len(panels) == 2:
      fig = plt.figure(figsize=(16, 7))
      fig.subplots_adjust(wspace=0.17, left=0.07, bottom=0.087, right=0.975, top=0.985)
  
  ax = {}
  for i_panel in range(len(panels)):
      
      panel = panels[i_panel]
      
      ax[panel] = fig.add_subplot(1,len(panels),i_panel+1)
      
      ts = t_BVI['t']
      
      #light curves
      if panel in ['lc', 'lc-zoom']:
          x_minmax, y_minmax = np.array([]), np.array([])
          for band in 'BVI':
              if panel == 'lc':
                  ax[panel].plot([ts, ts], [t_BVI[band]-t_BVI['e'+band], t_BVI[band]+t_BVI['e'+band]], '-', color=ecolor, lw=elw)
                  ax[panel].plot(ts, t_BVI[band], symbols[band], color='gray', markeredgecolor='k', ms=ms)
                  x_minmax = np.append(x_minmax, np.array([min(ts), max(ts)]))
                  y_minmax = np.append(y_minmax, np.array([min(t_BVI[band]-t_BVI['e'+band]), max(t_BVI[band]+t_BVI['e'+band])]))
              
              m = t_BVI[band][ts>=tmin]
              t = ts[ts>=tmin]
              m = m[t<=tmax]
              t = t[t<=tmax]
              
              if panel == 'lc-zoom':
                  x_minmax = np.append(x_minmax, np.array([min(t), max(t)]))
                  y_minmax = np.append(y_minmax, np.array([min(m), max(m)]))
              
              ax[panel].plot(t, m, symbols[band], color=colors[band], markeredgecolor='k', ms=ms, label='$'+band+'$')
      
      #Eh
      if panel in ['Eh', 'Eh-zoom']:
          
          if panel == 'Eh':
              x_minmax = t_EhBV_i.copy()
              ax[panel].plot([t_EhBV_i,t_EhBV_i], [EhBV_i-err_EhBV_i, EhBV_i+err_EhBV_i], '-', color=ecolor, lw=elw)
              ax[panel].plot(t_EhBV_i, EhBV_i, 'o', color='gray', markeredgecolor='k', ms=ms)
              y_minmax = [min(EhBV_i-err_EhBV_i),max(EhBV_i+err_EhBV_i)]
          
          EhBV_i     =     EhBV_i[t_EhBV_i>=tmin]
          err_EhBV_i = err_EhBV_i[t_EhBV_i>=tmin]
          t_EhBV_i   =   t_EhBV_i[t_EhBV_i>=tmin]
          EhBV_i     =     EhBV_i[t_EhBV_i<=tmax]
          err_EhBV_i = err_EhBV_i[t_EhBV_i<=tmax]
          t_EhBV_i   =   t_EhBV_i[t_EhBV_i<=tmax]
          ax[panel].plot([min(t_EhBV_i),max(t_EhBV_i)], [EhBV, EhBV], '-k')
          ax[panel].plot([min(t_EhBV_i),max(t_EhBV_i)], [EhBV-ssd_EhBV, EhBV-ssd_EhBV], '--k')
          ax[panel].plot([min(t_EhBV_i),max(t_EhBV_i)], [EhBV+ssd_EhBV, EhBV+ssd_EhBV], '--k')
          ax[panel].plot([t_EhBV_i,t_EhBV_i], [EhBV_i-err_EhBV_i, EhBV_i+err_EhBV_i], '-', color=ecolor, lw=elw)
          ax[panel].plot(t_EhBV_i, EhBV_i, 'o', color='lime', markeredgecolor='k', ms=ms)
          
          if panel == 'Eh-zoom':
              x_minmax = t_EhBV_i.copy()
              y_minmax = [min(EhBV_i-err_EhBV_i),max(EhBV_i+err_EhBV_i)]
          
          if i_panel == 0:  ax[panel].plot(-99, 99, ',w', label=sn)
          N_str = str(len(EhBV_i))
          N_text = '$N\!=\!$'+N_str
          ax[panel].plot(-99, 99, ',w', label=N_text)
          
          EhBV_str= str(EhBV)
          while len(EhBV_str)<5: EhBV_str = EhBV_str + '0'
          eEhBV_str= str(err_EhBV)
          while len(eEhBV_str)<5: eEhBV_str = eEhBV_str + '0'
          
          EhBV_text = '$E_{\mathrm{h}}(B\!-\!V)\!=\!$'+EhBV_str
          ax[panel].plot(-99, 99, ',w', label=EhBV_text)
          
          ssd_EhBV_str = str(round(ssd_EhBV,3))
          while len(ssd_EhBV_str)<5: ssd_EhBV_str = ssd_EhBV_str + '0'
          ssd_EhBV_text = '$\hat{\sigma}\!=\!$'+ssd_EhBV_str
          ax[panel].plot(-99, 99, ',w', label=ssd_EhBV_text)
          
          ax[panel].plot(-99, 99, ',w', label='$\sigma_{E_\mathrm{h}(B\!-\!V)}$='+eEhBV_str)
          
      #color-color plot
      if panel in ['c3', 'c3-zoom']:
          ts = np.atleast_1d(t_BVI['t'])
          Bs, eBs = np.atleast_1d(t_BVI['B']), np.atleast_1d(t_BVI['eB'])
          Vs, eVs = np.atleast_1d(t_BVI['V']), np.atleast_1d(t_BVI['eV'])
          Is, eIs = np.atleast_1d(t_BVI['I']), np.atleast_1d(t_BVI['eI'])
          x_sl, y_sl = np.array([]), np.array([])
          BVs, eBVs, VIs, eVIs = np.array([]), np.array([]), np.array([]), np.array([])
          for t, B, eB, V, eV, I, eI in zip(ts, Bs, eBs, Vs, eVs, Is, eIs):
          
              x, ex = B-V, np.sqrt(eB**2+eV**2)
              y, ey = V-I, np.sqrt(eV**2+eI**2)
              cov_xy = -eV**2
              
              BVs, eBVs, VIs, eVIs = np.append(BVs,x), np.append(eBVs,ex), np.append(VIs,y), np.append(eVIs,ey)
              
              if panel == 'c3':
                  plot_ellipse(ax[panel], x, ex, y, ey, cov_xy, ecolor, elw)
                  color = 'gray'
                  if t >= tmin and t <= tmax:
                      color='lime'
                      x_sl, y_sl = np.append(x_sl, x), np.append(y_sl, y)
                  ax[panel].plot(x, y, 'o', color=color, markeredgecolor='k', ms=ms)
                  
              if panel == 'c3-zoom':
                  if t >= tmin and t <= tmax:
                      color='lime'
                      x_sl, y_sl = np.append(x_sl, x), np.append(y_sl, y)
                      plot_ellipse(ax[panel], x, ex, y, ey, cov_xy, ecolor, elw)
                      ax[panel].plot(x, y, 'o', color=color, markeredgecolor='k', ms=ms)
          
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
          y_fit = 0.107 + 0.45*x_fit
          ax[panel].plot(x_fit, y_fit, ':b', zorder=0, label='Unreddened C3')
          
          
          if panel == 'c3':
              x_minmax = np.append(BVs-eBVs, BVs+eBVs)
              y_minmax = np.append(VIs-eVIs, VIs+eVIs)
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
          if i_panel == 0:  ax[panel].legend(loc='best', ncol=3, handletextpad=-0.3, columnspacing=0.1, edgecolor='none', title=sn)
          if i_panel == 1:  ax[panel].legend(loc='best', ncol=3, handletextpad=-0.3, columnspacing=0.1, edgecolor='none')
      if panel in ['Eh', 'Eh-zoom']:
          ax[panel].set_xlabel(time_label)
          ax[panel].set_ylabel('$E_{\mathrm{h}}(B\!-\!V)$')
          mpl.rcParams['legend.handlelength'] = 0.5
          ax[panel].legend(loc='best', ncol=1, handletextpad=-0.3, labelspacing=0.06, edgecolor='none')
      if panel in ['c3', 'c3-zoom']:
          ax[panel].set_xlabel('$B\!-\!V$')
          ax[panel].set_ylabel('$V\!-\!I$')
          mpl.rcParams['legend.handlelength'] = 1.5
          ax[panel].legend(loc='best', ncol=1, labelspacing=0.06, borderpad=0.2, handletextpad=0.0, edgecolor='none')
      
  if figure_name != '':  fig.savefig(figure_name+'.pdf', dpi=300)

