#generates a shell script for running multi variable parameter sweep
fid = open("driver_capped_final_tr1.sh","w")
import numpy as np
filexe = "./fil "
filini = "fil.ini"

folder = 'capped_pressure_model'

#s_tsyModelDir = tsydata
#s_tsyModel = tsymodel.dat

fid2 = open("./"+folder+"/data.midPointData.dat")

midData = []
for line in fid2.readlines():
	line = line.strip()
	columns = line.split()
	source = {}
	source['xe'] = columns[0]
	source['wb'] = columns[10]
	source['springkx'] = columns[9]
	midData.append(source)
fid2.close()
varNams = [
"f_apex",
"f_drive_omega_0",
"f_drive_omega_f",
"f_drive_coeff",
"f_drag_coeff",
"b_use_drag_force",
"b_use_driver_force",
"f_drive_tmin",
"f_drive_tmax",
"f_drag_tmin",
"f_drag_tmax",
"i_dimj",
"f_boundary",
"f_tMax",
"f_tInterv",
"f_pressureScale"
]

def myvalue(var,val):
	xe = float(val.get('xe'))
	ww = float(val.get('wb'))
	ts = 2.0*np.pi/ww
	if var == varNams[0]:
		my_value = val.get('xe')
	elif var == varNams[1]:
		my_value = str(0.5*ww)
	elif var == varNams[2]:
		my_value = str(2.0*ww)
	elif var == varNams[3]:
		my_value = str(1.0e-4*float(val.get('springkx')))
	elif var == varNams[4]:
		my_value = str(1.0e-4*float(val.get('springkx')))
	elif var == varNams[5]:
		my_value = str(0)
	elif var == varNams[6]:
		my_value = str(0)
	elif var == varNams[7]:
		my_value = str(0.0)
	elif var == varNams[8]:
		my_value = str(5.0*ts)
	elif var == varNams[9]:
		my_value = str(0.0)
	elif var == varNams[10]:
		my_value = str(5.0*ts)
	elif var == varNams[11]:
		my_value = str(501)
		#dimj
		#my_value = str(int(np.ceil(-500*(xe+2.0)/(12.0)+200)))
	elif var == varNams[12]:
		if abs(xe)<6.0:
			my_value = str(1.05)
		else:
			my_value = str(1.05+3.0*(abs(xe)-6.0)/8.0)
	elif var == varNams[13]:
		my_value = str(1000.0*ts)
	elif var == varNams[14]:
		my_value = str(0.1*ts)
	elif var == varNams[15]:
		my_value = str(1.0)
	return my_value

iii=0
for datum in reversed(midData):
	iii=iii+1
	line = filexe
	for var in varNams:
		line = line+"-i "+var+"="+myvalue(var,datum)+" "
	line = line+"-o "+"data_"+str(iii).zfill(3)+"_"+"apex"+"="+'%.3f'%(float(datum.get('xe')))
	line = line+" -i dontRunSim=0 -i plotPVG=0 -i buildTsyData=0"
	line = line+"\n"
	fid.write(line)
fid.close()