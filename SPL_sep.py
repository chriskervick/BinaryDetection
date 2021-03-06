#!/usr/bin/env python
from __future__ import absolute_import, unicode_literals, print_function
import math
import numpy as np
from scipy import integrate as integ
from scipy.special import erf
from scipy.special import hyp2f1
import sys
import time
import json
from numpy import pi, cos
from pymultinest.solve import solve
import os
from scipy.special import erfc 

#from mpi4py import MPI
Names = ['BootesI','BootesII','BootesIII','Coma','Crater','CanesI','CanesII','Draco','DracoII','Hercules','LeoI','LeoII','LeoIV','LeoV',\
    'PiscesII','SagII','SegueI','SegueII','SextansI','TriII','UrsaMajorI','UrsaMajorII','UrsaMinor','WillmanI','CetusII','ColumbaI',\
    'EridanusIII','Fornax','GrusI','GrusII','HoroI','HoroII','PhoenixII','PictorI','RetiII','RetiIII','TucanaII','TucanaIII','TucanaIV','TucanaV']


workdir,thisdir,pbf,pbg,spatialrundir,rundir,runnum = sys.argv[1:]

d = json.load(open(workdir+"/"+thisdir+"/dict_bf"+pbf+"bg"+pbg+".dat"))


dataloc = workdir+'/'+thisdir+'/rawxy'+runnum+'bf'+pbf+'bg'+pbg+'.dat'

rid = "run" + runnum +  "_bf" + pbf +"bg" + pbg

rundir = workdir + '/' + thisdir + '/' + rundir +'/'
spatialrundir = workdir + '/' + thisdir + '/' + spatialrundir +'/'

rp = d[Names[int(runnum)]]['Rp']
rmax = 10*rp


try: os.mkdir(rundir)
except OSError: pass
print(rundir)


                       
smin = d[Names[int(runnum)]]['lowlim']
s = np.loadtxt(rundir+"s_arr_"+rid+".dat")
s = s[s>smin]
smax = np.amax(s)
pi=math.pi;
numsamples=len(s)





def LogNormal(x,mu,sig,c,hc,erfscale):
    ans = c/(np.log(10)*x*np.sqrt(2*math.pi*sig**2)) * np.exp(-(np.log10(x)-mu)**2 / (2*sig**2))
    return(ans*0.5*(np.tanh((-erfscale*(x-hc))) + 1))


def myprior(cube):
    cube[0] = 10**(cube[0]*(10)-4)
    cube[1] = cube[1]*3.9-5
    return cube






def intbpl(smin,smax,g1,g2,sb,L): 
    hyp1 = hyp2f1((1+g1)*L,(g1-g2)*L,1+L+g1*L,-(smax/sb)**((1/L)))
    #print(hyp1)
    hyp2 = hyp2f1((1+g1)*L,(g1-g2)*L,1+L+g1*L,-(smin/sb)**((1/L)))
    #print(hyp2)
    c0 = (1/(1+g1))*(0.5**((-g1+g2)*L)) 
    c1 = smax*((smax/sb)**g1) 
    #print(c1)
    c2 = smin*(smin/sb)**g1
    #print(c2)
    ans = c0*(c1*hyp1 - c2*hyp2)
    return(ans)


def BPL(s,g1,g2,sb,L):
    ans = (s/sb)**g1 * (0.5*(1+(s/sb)**(1/L)))**((g2-g1)*L)
    return(ans)


def PlummerSep(s,a,Np):
    ata = (s*(2*a**2+s**2)*np.sqrt(4*a**2+s**2)) / (2*a**4 + 4*a**2*s**2 + s**4)
    ans = 4*Np**2*(-8*a**6*s+2*a**4*s**3+a**2*s**5+2*a**4*(a**2+s**2)*np.sqrt(4*a**2+s**2)*np.arctanh(ata)) / (4*a**2*s+s**3)**3
    return(ans*s)


def CrossTerm(s,a,rmax,Np,Nu):
    sqrtarg = (a**2+rmax**2)**2 + 2*(a-rmax)*(a+rmax)*s**2+s**4
    ans = (Np*Nu*(-a**2+rmax**2-s**2+np.sqrt(sqrtarg))) / (rmax**2*np.sqrt(sqrtarg))
    return(ans*s)


def tempUniCircle(s,rmax,Nu):
    s = s / (2*rmax)
    if s<1:
        ans = 16.0/math.pi * s * (np.arccos(s)-s*(1-s**2)**(0.5))
    else:
        ans = 0 
    return(Nu**2*ans /(2*rmax) )

def SPL(s,p):
    ans = s**(p)
    return(ans)
def SPLint(a,b,p):
    ans = (b**(p+1) - a**(p+1)) /(p+1)
    return(ans)

UniCircle = np.vectorize(tempUniCircle)


lowreslim = smin

#high = 2.0
spatialres = np.loadtxt(spatialrundir+'spatialparams'+rid+'.dat')
a = spatialres[0,0]
Np = spatialres[1,0]
Nu = spatialres[2,0]



print("SMIN IS " + str(smin) + " and SMAX IS " + str(smax))



#sb = 10**(3.97)*(4.84814*10**(-6))


def myloglike(cube):
    c,p=cube[0],cube[1]
    #c,g1,g2,sb,L = cube[0],cube[1],cube[2],cube[3],cube[4]
    ans = np.sum(np.log(PlummerSep(s,a,Np) + UniCircle(s,rmax,Nu)+2*CrossTerm(s,a,rmax,Np,Nu) + 0.5*SPL(s,p)*c))
    sepsinteg = integ.quad(PlummerSep,smin,smax,args=(a,Np))[0] + integ.quad(UniCircle,smin,smax,args=(rmax,Nu))[0] + 2*integ.quad(CrossTerm,smin,smax,args=(a,rmax,Np,Nu))[0]
    ans=ans-numsamples*np.log(0.5*c*SPLint(smin,smax,p) + sepsinteg)
    return(ans)

# number of dimensions our problem has
#parameters = ["c","g1","g2","sb",'L']
parameters = ["c","p"]
n_params = len(parameters)
# name of the output files

rid = rid+"SPL"
prefix = rundir + 'sepoutput' + rid
print(prefix)
# run MultiNest
result = solve(LogLikelihood=myloglike, Prior=myprior, n_dims=n_params, outputfiles_basename=prefix, verbose=True,resume=False)

print()
print('evidence: %(logZ).1f +- %(logZerr).1f' % result)
print()
print('parameter values:')
for name, col in zip(parameters, result['samples'].transpose()):
    print('%15s : %.3f +- %.3f' % (name, col.mean(), col.std()))


resforme = result['samples'].transpose()

c_val = np.mean(resforme[0,:])
c_std = np.std(resforme[0,:])
p_val =np.mean(resforme[1,:])
p_std =np.std(resforme[1,:])




resarray = np.array([[c_val,c_std],[p_val,p_std]])

np.savetxt(rundir+'sepparams'+rid+'.dat',resarray)






# make marginal plots by running:
# $ python multinest_marginals.py chains/3-
# For that, we need to store the parameter names:

#print("Took "+str(time.time()-start_time) + " seconds to get to JSON")
import json
with open('%sparams.json' % prefix, 'w') as f:
    json.dump(parameters, f, indent=2)
#print("Took "+str(time.time()-start_time) + " seconds to completely finish")

np.savetxt(rundir+'SEPTIME'+rid+'.dat',np.array([time.time()-start_time]))
