import numpy as np
import aplpy
import matplotlib.pyplot as plt
from astropy.io import fits
plt.clf()


targets = ['IRAC8uma','SOFIA7uma','SOFIA19uma','SOFIA31uma','SOFIA37uma','SOFIA37uma']
name = [r"$Spitzer \ 8 \mu m$",r"$SOFIA \  7 \mu m$",r"$SOFIA \ 19 \mu m$",r"$SOFIA \ 31 \mu m$",r"$SOFIA \ 37 \mu m$",r"$Herschel \ 70 \mu m$"]
beam = [2.0,3.0,2.5,2.8,2.9,5.1]
radio = [334.78949, 56.083436] # 37um peak
nstep = [7,5,7,7,7,7]  # number of contours
sigma = [1.06,3.20,3.67,9.45,8.22,8.22]  
bgn = [3, 4, 4, 4,5,5]  # sigma number of starting contour
label = ['(a)','(b)','(c)','(d)','(e)','(f)']
aper_fix = 3.84
aper_fix_mir1 = 4.6
aper_fix_mir3 = 4.6

fig = plt.figure(figsize=(10,10))
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)


i = 0

for y in [0.48,0.18]:
    for x in [0.09,0.39,0.69]:
        target = targets[i]
        data, hdr = fits.getdata(target+'.fits', 0, header = True)
        nanarr = np.where(np.isnan(data))
        data[nanarr] = np.median(data)
        figmin = np.min(data)
        figmax = np.max(data)
        data[nanarr] = figmax
        relative = data/figmax
        print(figmin/figmax)
        hdr['EQUINOX'] = 2000.0
        fits.writeto(target+'relative.fits', relative, hdr, clobber=True)
 

        
        data, hdr = fits.getdata(target+'relative.fits', 0, header = True)
        data0 = np.min(data)
        #data[np.where(data < sigma[i]*3./figmax)] = data0   # data below 3 sigma won't be shown
        fits.writeto('temp.fits', data, hdr, clobber=True)  

        if i == 5:
            f = aplpy.FITSFigure('temp.fits', subplot=[x-0.001,y+0.013,0.299,0.299],figure=fig)
        else:        
            f = aplpy.FITSFigure('temp.fits', subplot=[x,y,0.323,0.323],figure=fig)
            f.show_colorscale(cmap='CMRmap', vmin=1e-3,vmax=1,stretch='log')    


            bg0 = sigma[i]*bgn[i]    
            minl = np.log10(bg0)
            maxl = np.log10(figmax)
            dist = maxl-minl 
            step = 10**(dist/nstep[i])
            maxv = np.max(data)
            minv = maxv/step**nstep[i]
            n = np.arange(nstep[i])
            levels= step**n*minv
            print(levels)
            # llogf=np.append(background,llog2)
            f.show_contour(levels=levels,colors='k',linewidths=1.1)
            f.show_contour(levels=levels,colors='w',linewidths=0.8)


            # Label    
            fig.text(x+0.01,y+0.29,label[i],color='white')
            tex = r"${} \sigma = {:.2f} ; \ {:.2f}; \ {:.2f}$".format(bgn[i], bg0, step, figmax/1e3)
            fig.text(x+0.170,y+0.275,name[i],color='w',fontsize=12)  
            fig.text(x+0.135,y+0.025,tex,color='w',fontsize=9) 
            if target == 'SOFIA37uma':
                f.show_circles(radio[0], radio[1], aper_fix/3600., edgecolor='pink', facecolor='none', linestyle='dashed', linewidths=1.0)    # MIR2
                f.show_circles(334.78470,56.086256, aper_fix_mir1/3600., edgecolor='pink', facecolor='none', linestyle='dashed', linewidths=1.0)    # MIR1
                f.show_circles(334.78929,56.079328, aper_fix_mir3/3600., edgecolor='pink', facecolor='none', linestyle='dashed', linewidths=1.0)    # MIR3



            # Radio Source
            f.show_markers(radio[0],radio[1],edgecolor='k',facecolor='k', marker='+',s=80,zorder=5) # MIR2
            f.show_markers(334.78470,56.086256,edgecolor='k',facecolor='k', marker='+',s=80,zorder=5) # MIR1
            f.show_markers(334.78929,56.079328,edgecolor='k',facecolor='k', marker='+',s=80,zorder=5) # MIR3
            f.show_markers(334.7858333,56.08388889,edgecolor='w',facecolor='w', marker='+',linewidth=0.8,s=50,zorder=5)
            f.show_markers(334.7873917,56.08415833,edgecolor='w',facecolor='w', marker='+',linewidth=0.8,s=50,zorder=5) # 1.3mm MM2
            f.show_markers(334.7839792,56.08650278,edgecolor='w',facecolor='w', marker='+',linewidth=0.8,s=50,zorder=5) # MM1
            f.show_markers(334.7882417,56.08399167,edgecolor='w',facecolor='w', marker='+',linewidth=0.8,s=50,zorder=5) # MM3
            f.show_markers(334.7883292,56.08462778,edgecolor='w',facecolor='w', marker='+',linewidth=0.8,s=50,zorder=5) # MM4
            
            # Labels
            if i == 0:
                f.add_label(radio[0]+0.002,radio[1],'MIR2',color='white')
                f.add_label(334.78670,56.086256,'MIR1',color='white')
                f.add_label(334.79129,56.079328,'MIR3',color='white')
                f.add_label(334.7833333, 56.08388889,'3.4mm',color='white')
                f.add_label(334.7863917,56.08455833,'MM2',color='white')
                f.add_label(334.7819792,56.08650278,'MM1',color='white')
                f.add_label(334.7872417,56.08339167,'MM3',color='white')                
                f.add_label(334.7883292,56.08512778,'MM4',color='white')
            
           
            # Beam
            f.show_circles(f.pixel2world(100+1,100+1)[0],f.pixel2world(100+1,100+1)[1],beam[i]/7200.,edgecolor='none',facecolor='gray')
            # Outflow
            pa1 = 13./180.*np.pi
            if i == 0:            
                f.show_lines([np.array([[radio[0],(radio[0]+0.02*np.sin(pa1)/np.cos(radio[1]/180.*np.pi))],[radio[1],(radio[1]+0.02*np.cos(pa1))]])],edgecolor='k',linewidths=2.0)
                f.show_lines([np.array([[radio[0],(radio[0]+0.02*np.sin(pa1)/np.cos(radio[1]/180.*np.pi))],[radio[1],(radio[1]+0.02*np.cos(pa1))]])],edgecolor='c',linewidths=1.0)
                f.show_lines([np.array([[radio[0],(radio[0]+0.02*np.sin(pa1+np.pi)/np.cos(radio[1]/180.*np.pi))],[radio[1],(radio[1]+0.02*np.cos(pa1+np.pi))]])],linestyles='dashed',edgecolor='k',linewidths=2.0)
                f.show_lines([np.array([[radio[0],(radio[0]+0.02*np.sin(pa1+np.pi)/np.cos(radio[1]/180.*np.pi))],[radio[1],(radio[1]+0.02*np.cos(pa1+np.pi))]])],linestyles='dashed',edgecolor='r',linewidths=1.0)


            
            f.add_colorbar()
            f.colorbar.set_box([0.095,0.13,0.89,0.02],box_orientation='horizontal')
            #f.colorbar.set_location('top')
            #f.colorbar.set_width(0.1)
            f.colorbar.set_font(size='xx-small')
            f.colorbar.set_ticks([1.0,1e-1,1e-2,1e-3])
            f.colorbar.set_labels([1.0,1e-1,1e-2,1e-3])
        f.tick_labels.set_font(size='xx-small')
        #f.ticks.set_color('black')
        f.axis_labels.set_font(size='xx-small')
        f.tick_labels.set_xformat('hh:mm:ss')
        f.tick_labels.set_yformat('dd:mm:ss')  
        #f.ticks.set_xspacing(4.*15./3600.)
        f.tick_labels.hide()
        f.axis_labels.hide() 
        if i == 0 or i == 3:
            f.axis_labels.show_y()
            f.tick_labels.show_y()
        if i == 3 or i == 4 or i == 5:
            f.axis_labels.show_x()
            f.tick_labels.show_x()
        i+=1
#fig.text(0.42,0.845,'No SOFIA 7um or 11um data',color='k')   
fig.text(0.765,0.343,r"$No \ Herschel \ 70 \mu m \ data$",color='k')
fig.text(0.09,0.8,'IRAS22172')
fig.savefig('IRAS22172panel.eps',bbox_inches='tight')        
        
