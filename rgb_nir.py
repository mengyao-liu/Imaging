import numpy as np
import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits


plt.clf()
fig = plt.figure(figsize=(10,10))
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
target = 'IRAS22172'
imglist = ['UKIDSS_K_2.fits','UKIDSS_H_2.fits','UKIDSS_J_2.fits']
distance = {'S235':1800., 'IRAS22198':764., 'NGC2071':390., 'CepE':730., 'L1206':776., 'IRAS22172':2400., 'IRAS21391':750.}





#aplpy.make_rgb_cube(imglist,'thecube_NIR.fits')
#aplpy.make_rgb_image('thecube_NIR.fits',target+'_NIR_rgb.eps',stretch_r='arcsinh', stretch_g='arcsinh',stretch_b='arcsinh')
f = aplpy.FITSFigure('thecube_NIR_2d.fits', subplot=[0.1,0.1,0.8,0.8],figure=fig)
f.recenter(334.7869333, 56.08265556,0.0063999999)
f.show_rgb(target+'_NIR_rgb.eps')       
  

fsofia = fits.open('SOFIA37umarelative.fits')
#fsofia = fits.open('IRAC8umarelative.fits')
levs8 = [0.00668756,  0.01367551,  0.02796531,  0.05718679,  0.11694234,  0.23913759, 0.48901696]
levs37 = [0.09313977, 0.13073792, 0.18351348, 0.25759319, 0.36157699, 0.5075364, 0.71241589]
f.show_contour(fsofia,levels=levs37,colors='w',linewidths=0.8,zorder=4)
#f.show_contour(fsofia,levels=levs8,colors='w',linewidths=0.8,zorder=4)


f.show_markers(334.78949, 56.083436,edgecolor='k',facecolor='k', marker='+',s=400,linewidth=1.2,zorder=5) # MIR2
f.show_markers(334.78470,56.086256,edgecolor='k',facecolor='k', marker='+',s=400,linewidth=1.2,zorder=5) # MIR1
f.show_markers(334.78929,56.079328,edgecolor='k',facecolor='k', marker='+',s=400,linewidth=1.2,zorder=5) # MIR3
      
        
f.ticks.set_color('black')       
f.tick_labels.set_xformat('hh:mm:ss')
f.tick_labels.set_yformat('dd:mm:ss')
f.ticks.hide()
f.tick_labels.hide()
f.axis_labels.hide()

fig.text(0.45,0.91,target,color='k',fontsize=20)    


# show the scalebar
startpoint = f.pixel2world(263,32)
scalebar = np.array(([startpoint[0],startpoint[0]-10./3600.],[startpoint[1],startpoint[1]]))
f.show_lines([scalebar],edgecolor='w',linewidths=2.0)
fig.text(0.84,0.136,'10"',color='w',size='large')
length = 10./3600.*np.pi/180.*distance[target]
fig.text(0.825,0.114,'%0.2f'%length+'pc',color='w',size='large')
               


        
fig.savefig(target+'_NIR_37contour.eps',bbox_inches='tight')
#fig.savefig(target+'_NIR_8contour.eps',bbox_inches='tight')



