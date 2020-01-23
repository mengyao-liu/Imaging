import numpy as np
import aplpy
import matplotlib.pyplot as plt
import os
from astropy.io import fits


plt.clf()
fig = plt.figure(figsize=(30,30))
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
targets =['S235','IRAS22198','NGC2071','CepE','L1206','IRAS22172','IRAS21391']
imglist = ['SOFIA37umarelative.fits','SOFIA19umarelative.fits','IRAC8umarelative.fits']
distance = {'S235':1800., 'IRAS22198':764., 'NGC2071':390., 'CepE':730., 'L1206':776., 'IRAS22172':2400., 'IRAS21391':750.}
beam = [2.0,3.3,3.5]


i = 0
for y in [0.29,0.05]:  
    for x in [0.03,0.27,0.51,0.75]:
        if i == 7:
            break
        target = targets[i]
        #os.system('cd '+name+'/Rebinned')
        os.chdir('/Users/cici/Desktop/research/SOFIA/SOMA III/'+target)
        #aplpy.make_rgb_cube(imglist,'thecube.fits')
        aplpy.make_rgb_image('thecube.fits',target+'_rgb.eps',stretch_r='arcsinh', stretch_g='arcsinh',stretch_b='arcsinh',vmin_r=1e-2,vmin_g=1e-2,vmin_b=1e-2,vmax_r=1., vmax_g=1.,vmax_b=1.)   
        #aplpy.make_rgb_image('thecube.fits',target+'_rgb.eps',stretch_r='power',stretch_g='power',stretch_b='power',vmin_r=1e-2,vmin_g=1e-2,vmin_b=1e-2,vmax_r=1.0,vmax_g=1.0,vmax_b=1.0,exponent_r =0.5,exponent_g=0.5,exponent_b=0.5)   
        f = aplpy.FITSFigure('thecube_2d.fits', subplot=[x,y,0.21,0.21],figure=fig)
        f.show_rgb(target+'_rgb.eps')
        
        
        # show the scalebar
        startpoint = f.pixel2world(900,60)
        scalebar = np.array(([startpoint[0],startpoint[0]-10./3600.],[startpoint[1],startpoint[1]]))
        f.show_lines([scalebar],edgecolor='w',linewidths=2.0,zorder=10)
        fig.text(x+0.195,y+0.015,'10"',color='w',size='large')
        length = 10./3600.*np.pi/180.*distance[target]
        fig.text(x+0.19,y+0.006,'%0.2f'%length+'pc',color='w',size='large')
        
        # Beam
        f.show_circles(f.pixel2world(50,180)[0],f.pixel2world(50,180)[1],beam[0]/7200.,edgecolor='none',facecolor='blue')
        fig.text(x+0.02,y+0.035,r"$blue \ 8 \mu m$",color='w',size='large')
        f.show_circles(f.pixel2world(50,120)[0],f.pixel2world(50,120)[1],beam[1]/7200.,edgecolor='none',facecolor='green')
        fig.text(x+0.02,y+0.023,r"$green \ 19 \mu m$",color='w',size='large')
        f.show_circles(f.pixel2world(50,50)[0],f.pixel2world(50,50)[1],beam[2]/7200.,edgecolor='none',facecolor='red')
        fig.text(x+0.02,y+0.008,r"$red \ 37 \mu m$",color='w',size='large')
        
        
        #f.tick_labels.set_font(size='large')
        f.axis_labels.set_font(size='x-large')
        f.ticks.set_color('black')       
        f.tick_labels.set_xformat('hh:mm:ss')
        f.tick_labels.set_yformat('dd:mm:ss')
        fig.text(x+0.08,y+0.215,target,color='k',size='xx-large')
        f.axis_labels.hide()
        if i == 0 or i == 4:
            f.axis_labels.show_y()
        if i == 4 or i == 5 or i == 6 or i == 7:
            f.axis_labels.show_x()
            
        i+=1


fig.savefig('/Users/cici/Desktop/research/SOFIA/rgb/rgbsoma3_arcsinh.eps',bbox_inches='tight')
os.chdir('/Users/cici/Desktop/research/SOFIA/rgb')



