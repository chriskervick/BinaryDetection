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
pi=math.pi;
Names = ['BootesI','BootesII','BootesIII','Coma','Crater','CanesI','CanesII','Draco','DracoII','Hercules','LeoI','LeoII','LeoIV','LeoV',\
    'PiscesII','SagII','SegueI','SegueII','SextansI','TriII','UrsaMajorI','UrsaMajorII','UrsaMinor','WillmanI','CetusII','ColumbaI',\
    'EridanusIII','Fornax','GrusI','GrusII','HoroI','HoroII','PhoenixII','PictorI','RetiII','RetiIII','TucanaII','TucanaIII','TucanaIV','TucanaV']

workdir,thisdir,pbf,pbg,rundir,runnum = sys.argv[1:]


d = json.load(open(workdir+"/"+thisdir+"/dict_bf"+pbf+"bg"+pbg+".dat"))
dataloc = workdir+'/'+thisdir+'/rawxy'+runnum+'bf'+pbf+'bg'+pbg+'.dat'

rid = "run" + runnum +  "_bf" + pbf +"bg" + pbg

rundir = workdir + '/' + thisdir + '/' + rundir +'/'
try: os.mkdir(rundir)
except OSError: pass

# probability function, taken from the eggbox problem.
start_time = time.time()

bxy = np.loadtxt(dataloc)
bx = bxy[:,0]
by = bxy[:,1]
brs = np.sqrt(bx**2 + by**2)
rp = d[Names[int(runnum)]]['Rp']
brs = brs[brs<10*rp]
numSamples = len(brs)


#bans = CorrNum(bx,by,bins)






#r = np.loadtxt("10k_binary_rs_v2.dat")
Rmax = 10*rp


def Plummer(r,a,Np):
    ans = Np/(math.pi*a**2) * 1/(1+(r/a)**2)**2
    return(ans)
def Uniform(c,Rmax):
    return(1/(math.pi*Rmax**2) * c)

def myprior(cube):
    cube[0] = 10**(cube[0]*3)
    cube[1] = 10**(cube[1]*7.1-2)  
    cube[2] = 10**(cube[2]*7.1-2)  
    return cube

def myloglike(cube):
    a,Np,Nu = cube[0],cube[1],cube[2]
    ans = np.sum(np.log(Plummer(brs,a,Np) + Uniform(Nu,Rmax)))
    ans = ans - (Np*(1-a**2/(a**2+Rmax**2))) - Nu
    return(ans)


# number of dimensions our problem has
parameters = ["a", "Np","Nu"]
n_params = len(parameters)
# name of the output files
prefix = rundir+'output_' + rid

# run MultiNest
result = solve(LogLikelihood=myloglike, Prior=myprior, 
               n_dims=n_params, outputfiles_basename=prefix, verbose=True,resume=False)

print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
for name, col in zip(parameters, result['samples'].transpose()):
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))
    print(col.mean(),col.std())

resforme = result['samples'].transpose()

a_val = np.mean(resforme[0,:])
a_std = np.std(resforme[0,:])
Np_val =np.mean(resforme[1,:])
Np_std =np.std(resforme[1,:])
Nu_val = np.mean(resforme[2,:])
Nu_std = np.std(resforme[2,:])

resarray = np.array([[a_val,a_std],[Np_val,Np_std],[Nu_val,Nu_std]])
print(rundir+'spatialparams'+rid+'.dat')
np.savetxt(rundir+'spatialparams'+rid+'.dat',resarray)

# make marginal plots by running:
# $ python multinest_marginals.py chains/3-
# For that, we need to store the parameter names:

print("Took "+str(time.time()-start_time) + " seconds to get to JSON")
import json
with open('%sparams.json' % prefix, 'w') as f:
    json.dump(parameters, f, indent=2)
print("Took "+str(time.time()-start_time) + " seconds to completely finish")


np.savetxt(rundir+'SpatialTime'+rid+'.dat',np.array([time.time()-start_time]))


