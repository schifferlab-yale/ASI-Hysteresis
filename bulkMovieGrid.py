import glob
import numpy as np
import matplotlib as mpl


def decodeAttributes(string):
    dict={}
    desc=string.split(";")
    for attribute in desc:
        dict[attribute[0]] = float(attribute[1:])
    return dict["p"], dict["a"], dict["l"], dict["s"]




filePrefix="movieData/"
folders=sorted(glob.glob(filePrefix+"*.out"))
allSims=[]
for folder in folders:
    p,a,l,s=decodeAttributes(folder[len(filePrefix):-3])
    params={"p":p,"a":a,"l":l,"address":folder}
    allSims.append(params)



for spacing in np.array([280,300,320,340,360,380,400,420,440,512,1024])*1e-9:
    theseRuns=[sim for sim in allSims if abs(sim["a"]-spacing)<1e-9]
    print(theseRuns)

    break