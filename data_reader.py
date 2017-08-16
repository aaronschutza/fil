import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
from collections import defaultdict
import re
import sys

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
    for datum in fulldat:
        t = datum.get('t').tolist().pop()
        taxis.append(t)
        vx = datum.get('vx')
        vz = datum.get('vz')
        vxa = np.average(vx)
        vza = np.average(vz)
        vxavg.append(vxa)
        vzavg.append(vza)
        xdata = (datum.get('x')).tolist()
        zdata = (datum.get('z')).tolist()
        xmid.append(xdata[len(xdata)/2-1])
        zmid.append(zdata[len(zdata)/2-1])
        vxdata = (datum.get('vx')).tolist()
        vzdata = (datum.get('vz')).tolist()
        vxmid.append(vxdata[len(vxdata)/2-1])
        vzmid.append(vzdata[len(vzdata)/2-1])
    
    # N = len(taxis)
    # T = taxis[1]-taxis[0]
    # x = taxis
    # y = vxavg
    
    #y = vzavg
    #y = xmid-np.average(xmid)
    #y = zmid-np.average(zmid)
    #y = vxmid
    #y = vzmid

    N = len(taxis)/2
    T = taxis[1]-taxis[0]
    x = taxis[len(taxis)/2:]
    xavg = np.mean(xmid[len(taxis)/2:])
    #y = xmid[len(taxis)/2:]-xavg
    y = vxavg[len(taxis)/2:]

    yf = scipy.fftpack.fft(y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
    psyf = (2.0/N * np.abs(yf[:N//2])).tolist()
    
    #xf = xf[2:13]
    #psyf = psyf[2:13]
    
    spline = UnivariateSpline(xf, psyf-np.max(psyf)/2, s=0)
    roots = []
    max_value = max(psyf)
    max_index = psyf.index(max_value)
    omega_B = 2.0*np.pi*xf[max_index]
    ob.append(omega_B)
    for r in spline.roots():
        roots.append(2.0*np.pi*r)
    if len(roots)>2:
        print('more than one peak')
        switch1 = 0
        for i in range(0,len(roots)-1):
            if (roots[i+1]>omega_B) & (roots[i]<omega_B):
                root1.append(roots[i+1])
                root2.append(roots[i])
                switch1 = 1
        if switch1==0:
            root1.append(0.0)
            root2.append(0.0)
    elif len(roots)==2:
        root1.append(roots.pop())
        root2.append(roots.pop())
    else:
        root1.append(0.0)
        root2.append(0.0)
    fftx.append(psyf)
    fftt.append(xf.tolist())
    #xx.append(re.findall("\d+\.\d+",file))
    xx.append(str(abs(xavg)))
    f.close()

fid2 = open("./"+folder+"/data.midPointData.dat")

xxg = []
wbg = []
for line in fid2.readlines():
    line = line.strip()
    columns = line.split()
    xxg.append(abs(float(columns[0])))
    wbg.append(float(columns[10]))
fid2.close()


plt.plot(xx,ob,'bs',xx,root1,'rs',xx,root2,'gs',xxg,wbg,'mo')
plt.ylabel('Angular Frequency [1/s]')
plt.xlabel('Equatorial Crossing X [Re]')
plt.legend(['f_b TFC','FWHM upper','FWHM lower','f_b linear'])
plt.show()

f = open(folder+'_omegab.dat','w')


for i in range(0,len(xx)):
    f.write(xx[i]+' '+str(root1[i])+' '+str(ob[i])+' '+str(root2[i])+'\n')
f.close()

#quit()

cm = plt.cm.get_cmap('RdYlBu')
xxx=[];yyy=[];zzz=[];sss=[]
for i in range(0,len(xx)):
    for j in range(0,len(fftt[i])):
        xxx.append(xx[i])
        #yyy.append(fftt[i][j]/max(fftt[i]))
        yyy.append(2.0*np.pi*fftt[i][j])
        zzz.append(fftx[i][j]/max(fftx[i]))
        sss.append(np.floor(200*fftx[i][j]/max(fftx[i])))

plt.title('Power Spectrum')
sc = plt.scatter(xxx, yyy, c=zzz, cmap=cm,edgecolors='none',marker='o',alpha=0.5,s=sss)
clb = plt.colorbar(sc)
clb.set_label('Normalized Power')
plt.ylabel('Angular Frequency [1/s]')
plt.xlabel('Equatorial Crossing X [Re]')

plt.show()

