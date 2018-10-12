# Place holder for now. 
# Everything is copied from the original dash repo

def overall_offsets():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.24, left=0.12, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,4,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.17)
    fig.subplots_adjust(hspace=0.1)
    fs = 8

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    for j, file in enumerate(files):
        print file, '{}'.format(file.split('_')[1].split('.')[0])
        ax = fig.add_subplot(gs[j])
        if file.startswith('shifts_icxe15'):
            ax.set_xlim([-5.,70.])
            ax.set_ylim([-5.,20.])
        if file.startswith('shifts_icxe16'):
            ax.set_xlim([-5.,50.])
            ax.set_ylim([-10.,5.])
        if file.startswith('shifts_icxe17'):
            ax.set_xlim([-50.,5.])
            ax.set_ylim([-10.,5.])
        if file.startswith('shifts_icxe18'):
            ax.set_xlim([-40.,5.])
            ax.set_ylim([-5.,40.])
        ax.set_xlabel('$\Delta$ x [pix]', fontsize=fs, labelpad=0.1)
        if j == 0:
            ax.set_ylabel('$\Delta$ y [pix]', fontsize=fs)
        ax.set_title('{}'.format(file.split('_')[1].split('.')[0]), fontsize=fs)
            
        cc = 0.
        with open(file) as f:
            for line in f:
                if not line.startswith('#'):
                    data = line.split()
                    cc += 1.
                    color = scalarMap.to_rgba(cc)
                    #ax.arrow(0.,0., float(data[1]), float(data[2]), head_width=0., head_length=0., fc=color, ec = color)
                    ax.annotate("",xy=(float(data[1]), float(data[2])), xytext=(0,0), 
                        arrowprops=dict(arrowstyle='->', color=color))
        
    plt.show(block=False)
    fig.savefig('overall_offsets.png', dpi=200, transparent=False)


def overall_offsets_vs_distance():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.47, left=0.2, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.2)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=8)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    lines = ['-','--','-.',':']
    labels = ['COSMOS-15', 'COSMOS-16', 'COSMOS-17','COSMOS-18']
    
    for j, file in enumerate(files):
        print file
        root = '{}'.format(file.split('_')[1].split('.')[0])
        ax1.set_xlim([-0.5,9.])
        ax1.set_ylim([-0.005,0.15])
        ax1.set_xlabel('Commanded Offset Relative to First Position\n[arcmin]', fontsize=fs, labelpad=0.1)
        ax1.set_ylabel('Offset from Commanded Position\n[arcmin]', fontsize=fs)
        root = '{}'.format(file.split('_')[1].split('.')[0])
            
        data = table.read(file, format='ascii',  names=('file','x','y','rot','scale','x_rms','y_rms'))    
        small_off  = np.sqrt(data['x']**2 + data['y']**2)*0.12/60.
        #colors = (scalarMap.to_rgba(cc) for cc in range(len(small_off)))
        big_off = [0.0]
        drift_rate = [0.0]
        
        origin = fits.open(data['file'][0])
        x_orig, y_orig = origin[1].header['CRVAL1'], origin[1].header['CRVAL2'] 
        
        for k, file in enumerate(data['file'][1:]):
            flt = fits.open(file)
            x, y = flt[1].header['CRVAL1'], flt[1].header['CRVAL2']
            big_off.append(60.*np.sqrt((x_orig-x)**2 + (y_orig-y)**2))
            
            if k+1 < 4:
                t_exp = 255
            else:
                t_exp = 277
            drift = table.read(file.split('_')[0]+'_shifts.txt', format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
            drift_rate.append(25.*np.sqrt((drift['x'][0]-drift['x'][-1])**2 + (drift['y'][0]-drift['y'][-1])**2)/t_exp)
            
        print len(small_off), len(big_off)
        ax1.plot(big_off, small_off,linestyle=lines[j], color='0.5', label=labels[j], zorder=0)
        ax1.scatter(big_off, small_off, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        ax1.legend(loc='upper left', frameon=False, labelspacing=0.8, fontsize=9)        
        
        ax2.set_xlim([-0.5, 9.])
        ax2.set_ylim([-0.01,0.5])
        ax2.set_xlabel('Commanded Offset Relative to First Position\n[arcmin]', fontsize=fs, labelpad=0.1)
        ax2.set_ylabel('Drift Rate During Observations\n[pix per 25 seconds]', fontsize=fs)
        
        ax2.plot(big_off, drift_rate,linestyle=lines[j], color='0.5', zorder=0)
        ax2.scatter(big_off, drift_rate, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        
        
        
    plt.show(block=False)
    fig.savefig('overall_offsets_vs_distance.png', dpi=200, transparent=False)


def overall_offsets_vs_time():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    files = glob.glob('shifts_icxe*010.txt')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.47, left=0.2, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,2,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.2)
    fig.subplots_adjust(hspace=0.2)
    fs = 10
    matplotlib.rc('xtick',labelsize=fs)
    matplotlib.rc('ytick',labelsize=fs)

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=8)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    time = np.array([0, 368, 713, 1060, 1405, 1777, 2147, 2519]) + 333.
    lines = ['-','--','-.',':']
    labels = ['COSMOS-15', 'COSMOS-16', 'COSMOS-17','COSMOS-18']
    
    for j, file in enumerate(files):
        root = '{}'.format(file.split('_')[1].split('.')[0])
        ax1.set_xlim([-5., 3100.])
        ax1.set_ylim([-0.005,0.15*60.])
        ax1.set_xlabel('Time Since Beginning of Orbit\n[seconds]', fontsize=fs)
        ax1.set_ylabel('Offset from Commanded Position\n[arcsec]', fontsize=fs)
        root = '{}'.format(file.split('_')[1].split('.')[0])
            
        data = table.read(file, format='ascii',  names=('file','x','y','rot','scale','x_rms','y_rms'))    
        small_off  = np.sqrt(data['x']**2 + data['y']**2)*0.12
        #colors = (scalarMap.to_rgba(cc) for cc in range(len(small_off)))
        drift_rate = [0.0]
                
        for k, file in enumerate(data['file'][1:]):
            if k+1 < 4:
                t_exp = 255
            else:
                t_exp = 277
            drift = table.read(file.split('_')[0]+'_shifts.txt', format='ascii', 
                names=('file','x','y','rot','scale','x_rms','y_rms'))
            drift_rate.append(25.*np.sqrt((drift['x'][0]-drift['x'][-1])**2 + (drift['y'][0]-drift['y'][-1])**2)/t_exp)
            
        ax1.plot(time, small_off, linestyle=lines[j], color='0.5', label=labels[j], zorder=0)
        ax1.scatter(time, small_off, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        ax1.legend(loc='upper left', frameon=False, labelspacing=0.8, fontsize=9)
        
        
        ax2.set_xlim([-5., 3100.])
        ax2.set_ylim([-0.01,0.5])
        ax2.set_xlabel('Time Since Beginning of Orbit\n[seconds]', fontsize=fs)
        ax2.set_ylabel('Drift Rate During Pointing\n[pix per 25 seconds]', fontsize=fs)
        
        ax2.plot(time, drift_rate,linestyle=lines[j], color='0.5', zorder=0)
        ax2.scatter(time, drift_rate, c=range(len(small_off)), cmap='jet', s=35., edgecolors='black', alpha=0.7)
        
    plt.show(block=False)
    fig.savefig('overall_offsets_vs_time.png', dpi=200, transparent=False)
    fig.savefig('overall_offsets_vs_time.pdf', dpi=200, transparent=False)

        
def gyro_drift():
    
    import glob
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import threedhst
    
    files = glob.glob('icxe*010_asn.fits')
    
    fig = unicorn.catalogs.plot_init(xs=10.,aspect=0.24, left=0.12, right=0.05, bottom=0.15, top=0.15, NO_GUI=False)    
    gs = gridspec.GridSpec(1,4,top=0.9, bottom=0.15)
    fig.subplots_adjust(wspace=0.17)
    fig.subplots_adjust(hspace=0.1)
    fs = 8

    jet = cm = plt.get_cmap('jet_r')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
    
    for j, file in enumerate(files):
        print file
        ax = fig.add_subplot(gs[j])
        ax.set_xlabel('$\Delta$ x [pix]', fontsize=fs, labelpad=0.1)
        if j == 0:
            ax.set_ylabel('$\Delta$ y [pix]', fontsize=fs)
        ax.set_title('{}'.format(file.split('_')[0]), fontsize=fs)
        ax.set_xlim([-3., 3.])
        ax.set_ylim([-3., 3.])
        
        asn = threedhst.utils.ASNFile(file)
        cc = 1.
        for exp in asn.exposures[1:]:
            
            data = table.read(exp+'_shifts.txt', format='ascii', names=('file','x','y','rot','scale','x_rms','y_rms'))
            cc += 1.
            color = scalarMap.to_rgba(cc)
            ax.plot(data['x'], data['y'], '-', color=color, alpha=0.5, linewidth=1.5)
            ax.plot(data['x'], data['y'], 'o', color=color, markersize=3., alpha=0.5, markeredgecolor=None)
            
    plt.show(block=False)
    fig.savefig('gyro_drift.png', dpi=200, transparent=False)                    

def footprints_plot(root='icxe15010'):
    
    import unicorn.survey_paper as sup
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    
    if root == 'icxe15010':
        aspect = 1.75
        xlim = [150.265, 150.157]
        ylim = [2.45, 2.64]
        xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$']
        xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True)]
        yticklab = [r'$+02^\circ30^\prime00^{\prime\prime}$',r'$+02^\circ35^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 30, 00, hours=False),sup.degrees(2, 35, 00, hours=False)]
        label = 'COSMOS-15'
        factor=10.

    if root == 'icxe16010':
        aspect=0.9
        xlim = [150.265, 150.1]
        ylim = [2.607, 2.74]
        xticklab = [r'$10^\mathrm{h}01^\mathrm{m}00^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$']
        xtickv = [sup.degrees(10,01,00, hours=True),sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True)]
        yticklab = [r'$+02^\circ38^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$', r'$+02^\circ42^\prime00^{\prime\prime}$', r'$+02^\circ44^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 38, 00, hours=False),sup.degrees(2, 40, 00, hours=False),sup.degrees(2, 42, 00, hours=False),sup.degrees(2, 44, 00, hours=False)]
        label='COSMOS-16'
        factor=20.
    
    if root == 'icxe17010':
        aspect=1.4
        xlim = [150.2, 150.06]
        ylim = [2.52, 2.72]
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}45^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}15^\mathrm{s}$']
        xtickv = [sup.degrees(10,00,45, hours=True),sup.degrees(10,00,30, hours=True),sup.degrees(10,00,15, hours=True)]
        yticklab = [r'$+02^\circ35^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False)]
        label='COSMOS-17'
        factor=240.

    if root == 'icxe18010':
        aspect=1.577
        xlim = [150.14, 150.01]
        ylim = [2.53, 2.735]
        xticklab = [r'$10^\mathrm{h}00^\mathrm{m}30^\mathrm{s}$', r'$10^\mathrm{h}00^\mathrm{m}20^\mathrm{s}$',r'$10^\mathrm{h}00^\mathrm{m}10^\mathrm{s}$']
        xtickv = [sup.degrees(10,00,30, hours=True),sup.degrees(10,00,20, hours=True),sup.degrees(10,00,10, hours=True)]
        yticklab = [r'$+02^\circ35^\prime00^{\prime\prime}$',r'$+02^\circ40^\prime00^{\prime\prime}$']
        ytickv = [sup.degrees(2, 35, 00, hours=False),sup.degrees(2, 40, 00, hours=False)]
        label='COSMOS-18'
        factor=240.
    
    
    
    fig = unicorn.catalogs.plot_init(square=True, xs=5., aspect=aspect, 
        fontsize=8, left=0.18, right=0.02, bottom=0.10, top=0.10)
    ax = fig.add_subplot(111)
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=9)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)    
    
    reg_file = root+'_asn.reg'
    
    poly = []
    with open(reg_file) as f:
        for line in f:
            if not line.startswith('fk5'):
                region = line.split('#')[0]
                poly.append(sup.polysplit(region=region, get_shapely=True))

    shifts = table.read('shifts_{}.txt'.format(root), format='ascii', 
        names=('file','x','y','rot','scale','x_rms','y_rms'))
        
    cc = 0
    xcen_all = []
    ycen_all = []
    for j,(pp, x_off, y_off, file) in enumerate(zip(poly, shifts['x'], shifts['y'], shifts['file'])):
        cc += 1.
        color = scalarMap.to_rgba(cc)
        x, y = pp.exterior.xy
        flt = fits.open(file)
        xcen = flt[1].header['CRVAL1O']
        ycen = flt[1].header['CRVAL2O']
        x_off = (flt[1].header['CRVAL1B']-flt[1].header['CRVAL1O'])*20.
        y_off = (flt[1].header['CRVAL2B']-flt[1].header['CRVAL2O'])*20.
        print file, xcen, xcen+x_off, ycen, ycen+y_off
        #xcen = (np.mean(x[:-1]))
        #ycen = (np.mean(y[:-1]))
        xcen_all.append(xcen)
        ycen_all.append(ycen)
        ax.plot(x,y,'-', color=color)
        #ax.annotate("",xy=(xcen+(x_off*0.12)/factor, ycen+(y_off*0.12)/factor), xytext=(xcen, ycen), 
        #    arrowprops=dict(arrowstyle='->', color=color))
        #ax.plot([xcen, xcen+x_off], [ycen, ycen+y_off], '-')
        ax.annotate("",xy=(xcen+x_off, ycen+y_off), xytext=(xcen, ycen), 
            arrowprops=dict(arrowstyle='->', color=color))

    ax.plot(xcen_all, ycen_all, '+:', markersize=10., color='0.5', alpha=0.5) 
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)     
    ax.set_xticklabels(xticklab)
    xtick = ax.set_xticks(xtickv)
    ax.set_yticklabels(yticklab)
    ytick = ax.set_yticks(ytickv)
    ax.set_title(label)       
    plt.show(block=False)
    
    fig.savefig('footprint_{}.png'.format(label.lower()), dpi=200, transparent=False)          
