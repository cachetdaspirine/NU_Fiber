import numpy as np
from numpy.random import default_rng
import math
import sys
import os

#sys.path.append('/home/hugo/Extra_Module_py')
sys.path.append('/home/hleroy/Simulation/Extra_Module_py')
from Numeric_Fiber_Energy import *
from Numeric_Hex_Energy import *
from Functions import *
import Conversion as Conv
import RandomParticleObject as RPO
import RandomParticleFunctions_v4 as RPF
import RandSyst as RSys
import Shape as Sh
import MeasurePoisson as MP



Nmax = 1000
Wmax = 10

def Get_GammaMax(bf,bh,P):
    dGamma = 0.01
    Neq=0
    while Neq < bh.Nmax:
        P.reComputeJ(P.Gamma+dGamma)
        Neq,E = bh.Get_Best_Disk(P)
    GammaHex = P.Gamma
    P.reComputeJ(0.)
    Weq=0
    while Weq < bf.Wmax-1:
        P.reComputeJ(P.Gamma+dGamma)
        Weq,E = bf.Get_Best_Fiber(P,type=1)
    return max(GammaHex,P.Gamma)

# Import the data
SimNum = int(sys.argv[1])
#Start = 4000
Nparticles = 10
Data = np.loadtxt('Matrix.data',dtype=float)
Data = Data[Data[:,1].argsort()]
dN = Data.shape[0]//100
rng = default_rng()
indexs = rng.choice(100,size=10,replace=False)
#Data = Data[Data.shape[0]-(SimNum+1)*Nparticles-Start:Data.shape[0]-(SimNum)*Nparticles-Start]
Data = Data[dN*SimNum:dN*(SimNum+1)][indexs]
Matrices = np.zeros((Nparticles,6),dtype=float)
Matrices[:,:-1] = Data

for n,ligne in enumerate(Matrices):
    Mc,rho0,e1,e2,seed = RPF.RandomParticle(int(ligne[-2]))
    Matrices[n,-1] = MP.ComputePoissonRatio(Mc,rho0)


# Generate the data
BestDiskLL=list()
BestFiberLL1=list()
Parameters = list()
seeds = np.zeros(Matrices.shape[0],dtype=np.int_)
DE = list()
for num, Matrice in enumerate(Matrices):
    BestDiskL=list()
    BestFiberL1 = list()
    Mc,rho0,e1,e2,seed = RPF.RandomParticle(seed=int(Matrice[-2]))
    Param = Conv.MatrixToContinuum(Mc,rho0,e1,e2,Gamma=0.)
    bh = BD(Nmax,Param,Expansion=True,Mc=Mc,q0=rho0)
    bf = BF(Wmax,Param,Expansion=True,Mc=Mc,q0=rho0)
    nu = Matrice[-1]#MP.ComputePoissonRatio(Mc,rho0)#ComputePoissonRatio(Particle.Mc,Particle.rho0)
    Ell = Matrice[1]#MeasureL(Mc,rho0)#Particle.length
    GammaMax = Get_GammaMax(bf,bh,Param)
    for gamma in np.linspace(0.,GammaMax,100):
        Param = Conv.MatrixToContinuum(Mc,rho0,e1,e2,Gamma=gamma)
        BestDiskL.append(bh.Get_Best_Disk(Param))
        BestFiberL1.append(bf.Get_Best_Fiber(Param,type=1))
    BestDiskL = np.array(BestDiskL)
    BestFiberL1 = np.array(BestFiberL1)
    DE.append(BestDiskL[:,1]-BestFiberL1[:,1])
    BestDiskLL.append(BestDiskL)
    BestFiberLL1.append(BestFiberL1)
    Parameters.append([nu,Ell,Matrice[2],Matrice[3],GammaMax])
    seeds[num] = int(Matrice[-2])

DE = np.array(DE)
BestFiberLL1 = np.array(BestFiberLL1)
BestDiskLL = np.array(BestDiskLL)
Parameters=np.array(Parameters)



#sort by Poisson ratio
BestDiskLL = BestDiskLL[Parameters[:,0].argsort()]
BestFiberLL1 = BestFiberLL1[Parameters[:,0].argsort()]
DE=DE[Parameters[:,0].argsort()]
Parameters = Parameters[Parameters[:,0].argsort()]

np.save('Disks_'+str(SimNum),BestDiskLL,allow_pickle=True)
np.save('Fiber_'+str(SimNum),BestFiberLL1,allow_pickle=True)

Aggregates=list()
for n in range(BestFiberLL1.__len__()):
    BFL1 = BestFiberLL1[n]#BFSorted1[n][0]
    BDL = BestDiskLL[n]
    BestAgg = np.zeros((BDL.shape[0],2),dtype=float)
    for i in range(BDL.shape[0]):
        if BDL[i,1]<BFL1[i,1] :#and BDL[i,1]<BFL2[i,1]:
            BestAgg[i] = [0,BDL[i,0]]
        else :
            BestAgg[i] = [BFL1[i,0],0]
    Aggregates.append(BestAgg)
np.save('Seeds_'+str(SimNum),seeds,allow_pickle=True)
np.save('DeltaE_'+str(SimNum),DE,allow_pickle=True)
np.save('Aggregates_'+str(SimNum),Aggregates,allow_pickle=True)
np.save('Parameters_'+str(SimNum),Parameters,allow_pickle=True)
