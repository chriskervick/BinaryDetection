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
brs = brs[brs<10*rp]






def f(s,a,I0):
    ans = 1/(s**3*(4*a**2 + s**2)**(5/2)) *4*a**6*I0**2 *math.pi**2*(s*(-2*a**2+s**2)*np.sqrt(4*a**2+s**2) - 2*a**2*(a**2+s**2)*(2*np.log(a)+np.log(-s + np.sqrt(4*a**2+s**2)) - np.log(s**2*(s + np.sqrt(4*a**2+s**2)) + a**2*(3*s+np.sqrt(4*a**2+s**2)))))
    return(s*ans)




slenmax = 200000

if(len(bx)**2 < 200):
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
    seps = 10000*np.ones(slenmax)
    maxval = np.amax(seps)
    maxidx = np.argmax(seps)
    for i in range(0,len(bx)):
        for j in range(i+1,len(bx)):
            thisdist = np.sqrt((bx[i]-bx[j])**2 + (by[i]-by[j])**2)
            if(thisdist<maxval):
                seps[maxidx]=thisdist
                maxval = np.amax(seps)
                maxidx = np.argmax(seps)
    ts = seps[seps<10000]
    s=ts[ts>0]
    print("SECOND LOOP, len(s) = " + str(len(s)))
    print("The highest separation we have reached is " + str(np.amax(s)))
    print("Just make sure it's more than 0.1pc")

np.savetxt(rundir+"s_arr_"+rid+".dat",s)
print("Took "+str(time.time()-start_time) + " seconds to completely finish")
np.savetxt(rundir+"MS_TIME"+rid+".dat",np.array([time.time()-start_time]))

