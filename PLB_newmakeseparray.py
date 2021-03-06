#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
import math
import numpy as np
from scipy import integrate as integ
import sys
import time
from numpy import pi, cos
from pymultinest.solve import solve
import os
from scipy.spatial.distance import cdist as cdist 
import json
#from mpi4py import MPI

Names = ['BootesI','BootesII','BootesIII','Coma','Crater','CanesI','CanesII','Draco','DracoII','Hercules','LeoI','LeoII','LeoIV','LeoV',\
    'PiscesII','SagII','SegueI','SegueII','SextansI','TriII','UrsaMajorI','UrsaMajorII','UrsaMinor','WillmanI','CetusII','ColumbaI',\
    'EridanusIII','Fornax','GrusI','GrusII','HoroI','HoroII','PhoenixII','PictorI','RetiII','RetiIII','TucanaII','TucanaIII','TucanaIV','TucanaV']




workdir,thisdir,pbf,pbg,spatialrundir,rundir,runnum = sys.argv[1:]



d = json.load(open(workdir+"/"+thisdir+"/dict_bf"+pbf+"bg"+pbg+".dat"))
dataloc = workdir+'/'+thisdir+'/rawxy'+runnum+'bf'+pbf+'bg'+pbg+'.dat'

rid = "run" + runnum +  "_bf" + pbf +'bg'+pbg

rundir = workdir + '/' + thisdir + '/' + rundir +'/'
spatialrundir = workdir + '/' + thisdir + '/' + spatialrundir +'/'





# probability function, taken from the eggbox problem.
start_time = time.time()



print(rundir)
try: os.mkdir(rundir)
except OSError: pass


#prevparams = np.loadtxt(spatialrundir+'spatialparams'+rid+'.dat')
#pa = prevparams[0,0]
#pI0 = prevparams[1,0]


bxy = np.loadtxt(dataloc)
bx = bxy[:,0]
by = bxy[:,1]
brs = np.sqrt(bx**2 + by**2)
rp = d[Names[int(runnum)]]['Rp']
bx = bx[brs<10*rp]
by = by[brs<10*rp]
bxy = bxy[brs<10*rp]
SEPLIM = 2.0

def RedoDiff(y):
    count = 0
    (s1,s2) = np.shape(y)
    t = np.zeros(s1*s2)
    for i in range(0,s1):
        for j in range(0,s2):
            t[count] = y[i,j]
            count+=1
    tt = t[t>0]
    ttt = tt[tt<SEPLIM]
    return(ttt)

def RedoSim(y):
    count = 0
    (s1,s2) = np.shape(y)
    t = np.zeros(s1*s2)
    for i in range(0,s1):
        for j in range(i+1,s2):
            t[count] = y[i,j]
            count+=1
    tt = t[t>0]
    ttt = tt[tt<SEPLIM]
    return(ttt)

cutsize = 500
sn = int(np.floor(len(bxy)/cutsize))

if(len(bxy) < cutsize):
    tbraws = np.zeros(len(bx)**2)
    count = 0
    for i in range(0,len(bx)):
        if i%100 == 0:
            print(i)
        for j in range(i+1,len(by)):
            tbraws[count] = np.sqrt((bx[i]-bx[j])**2 + (by[i]-by[j])**2)
            count+=1
    s = tbraws[tbraws>0]
    print("FIRST LOOP, len(s) = " + str(len(s)))
else:
    splitdat = np.array_split(bxy,sn)
    s = np.zeros(0)
    for i in range(0,sn):
        for j in range(i,sn):
            if i==j:
                y = cdist(splitdat[i],splitdat[i])
                tt = RedoSim(y)
                s = np.append(s,tt)
                del y 
                del tt
            else:
                y = cdist(splitdat[i],splitdat[j])
                tt = RedoDiff(y)
                s = np.append(s,tt)
                del y
                del tt







np.savetxt(rundir+"s_arr_"+rid+".dat",s)
print("Took "+str(time.time()-start_time) + " seconds to completely finish")
np.savetxt(rundir+"MS_TIME"+rid+".dat",np.array([time.time()-start_time]))

