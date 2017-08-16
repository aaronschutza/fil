import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from collections import defaultdict
import re
import sys
import matplotlib.animation as animation

def getdatafromfile():
    dimi=np.fromfile(f,dtype='int32',count=1)
    if not dimi:
        return [],[],1
    else:
        dimj=np.fromfile(f,dtype='int32',count=1)
        junk=np.fromfile(f,dtype='int32',count=1)
        time=np.fromfile(f,dtype='float64',count=1)
        step=np.fromfile(f,dtype='int32',count=1)
        data=np.reshape(np.fromfile(f,dtype='float64',count=dimi*dimj),(dimi,dimj),order='F')
        return time,data,0

folder = 'capped_pressure_model/'
filenames = []
import os
for file in os.listdir(folder):
    if file.endswith(".bin"):
        filenames.append(file)
filenames.sort()

ob = []
root1 = []
root2 = []
fftx = []
fftt = []
xx = []
from scipy.interpolate import UnivariateSpline
def make_norm_dist(x, mean, sd):
    return 1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x - mean)**2/(2*sd**2))

for file in filenames:
    print 'Working on '+file
    fulldat = []
    f = open(folder+file, 'rb')
    eof = 0
    while (eof<1):
        time,dat,eof=getdatafromfile()
        if (eof<1):
            source = defaultdict(list)
            #x = 1, z = 2, vx = 3, vz = 4, d = 5, p = 6, t = 7, b = 8, mf = 9, pdg = 10
            source['x'] = dat[0]
            source['z'] = dat[1]
            source['vx'] = dat[2]
            source['vz'] = dat[3]
            source['d'] = dat[4]
            source['p'] = dat[5]
            source['t'] = dat[6]
            source['b'] = dat[7]
            source['mf'] = dat[8]
            source['pdg'] = dat[9]
            source['t'] = time
            fulldat.append(source)
    taxis = []
    vxavg = []
    vzavg = []
    xmid = []
    zmid = []
    vxmid = []
    vzmid = []
    iii=0
    write_to_text=1
    if write_to_text==1:
        ftext = open(folder+file[0:-4]+'.txt','w')
        print 'writing text to '+folder+file[0:-4]+'.txt'
    for datum in fulldat:
        t = datum.get('t').tolist().pop()
        taxis.append(t)
        xdata = (datum.get('x'))
        dimj = len(xdata)
        zdata = (datum.get('z'))
        if write_to_text==1:
            ftext.write(str(dimj)+'\n')
            ftext.write(str(t)+'\n')
            for i in range(0,dimj):
                ftext.write(str(xdata[i])+' '+str(zdata[i])+'\n')
        if iii == 0:
            allxdata = xdata
            allzdata = zdata
            iii=1
        else:
            allxdata = np.vstack((allxdata, xdata))
            allzdata = np.vstack((allzdata, zdata))
    fig = plt.figure()
    ax = plt.axes(xlim=(-3.0, 0.0), ylim=(-2.0, 2.0))
    scat = ax.scatter([],[])
    def animate(i):
        plt.title('T='+str(taxis[i]))
        dddata = np.vstack((allxdata[i], allzdata[i])).transpose()
        scat.set_offsets(dddata)
        return scat,
    def init():
        scat.set_offsets([])
        return scat,
    plt.xlabel('X [Re]')
    plt.ylabel('Z [Re]')
    line_ani = animation.FuncAnimation(fig, animate, init_func=init, frames=range(len(taxis)),interval=50, blit=True)
    line_ani.save(folder+file+'.mp4')

    if write_to_text==1:
        ftext.close()

    # Y = np.asarray(taxis)
    # X = np.linspace(0, 1, dimj)
    # Z = vvt
    # Z = Z.reshape((len(Y),len(X)))
    # levels = np.linspace(0, 1, 40)
    # fig = plt.contourf(X, Y, Z)
    # plt.xlabel('Normalized position along filament')
    # plt.ylabel('Time [s]')
    # cbar = plt.colorbar(ticks=[0,1])
    # plt.suptitle(file)
    # plt.show()
