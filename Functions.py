import numpy as np
import sys
import os
import copy
sys.path.append('/home/hleroy/Simulation/Extra_Module_py')
#sys.path.append('/home/hugo/Extra_Module_py')
import RandSyst as RSys
import Shape as Sh

Name = str(sys.argv[1])

def DistanceFromEdge(Array):
    Res = dict()
    box = copy.copy(Array)
    A,B = np.where(Array==1)
    Index = set(list(zip(A,B)))
    ind = 0
    while Index.__len__()!=0:
        toremove=set()
        for ij in Index:
            if Sh.Get_Neighbors(box,ij,Free=True,ParticleType='Hexagon').__len__()!=0:
                Res[ij] = ind
                toremove.add(ij)
        #print(toremove)
        #input()
        for ij in toremove:
            box[ij]=0
        Index-=toremove
        ind+=1
    return Res
def MeasureL(Mc,rho0):
    Array = Sh.Parallel(18,ParticleType='Hexagon')
    S = RSys.System(Mc,rho0,Array)
    S.PrintPerSite(Name = Name,Extended=True)
    EnergyData = np.loadtxt(Name)
    Ranking = DistanceFromEdge(Array)
    os.system('rm '+Name)
    Res = np.zeros((max(Ranking.values())+1,2),dtype=float)
    for ligne in EnergyData:
        Res[Ranking[(int(ligne[-2]),int(ligne[-1]))],0] += ligne[-3]
        Res[Ranking[(int(ligne[-2]),int(ligne[-1]))],1] += 1
    #Res = Res[1:]
    Res[:,0] = Res[:,0]/Res[:,1]
    return sum((max(Res[:,0])-Res[:,0])/(max(Res[:,0])-min(Res[:,0])))
