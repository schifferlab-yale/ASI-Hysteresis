import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def read_mumax3_table(filename):
    """Puts the mumax3 output table in a pandas dataframe"""

    from pandas import read_table
    
    table = read_table(filename)
    table.columns = ' '.join(table.columns).split()[1::2]
    
    return table


def read_mumax3_ovffiles(outputdir):
    """Load all ovffiles in outputdir into a dictionary of numpy arrays 
    with the ovffilename (without extension) as key"""
    
    from subprocess import run, PIPE, STDOUT
    from glob import glob
    from os import path
    from numpy import load

    # convert all ovf files in the output directory to numpy files
    p = run(["mumax3-convert","-numpy",outputdir+"/*.ovf"], stdout=PIPE, stderr=STDOUT)
    if p.returncode != 0:
        print(p.stdout.decode('UTF-8'))

    # read the numpy files (the converted ovf files)
    fields = {}
    for npyfile in glob(outputdir+"/*.npy"):
        key = path.splitext(path.basename(npyfile))[0]
        fields[key] = load(npyfile)
    
    return fields

def run_mumax3(script, name, verbose=False):
    """ Executes a mumax3 script and convert ovf files to numpy files
    
    Parameters
    ----------
      script:  string containing the mumax3 input script
      name:    name of the simulation (this will be the name of the script and output dir)
      verbose: print stdout of mumax3 when it is finished
    """
    
    from subprocess import run, PIPE, STDOUT
    from os import path
    import os

    scriptfile = name + ".txt" 
    outputdir  = name + ".out"
    

    # write the input script in scriptfile
    with open(scriptfile, 'w' ) as f:
        f.write(script)
    
    # call mumax3 to execute this script
    p = run(["mumax3","-f",scriptfile], stdout=PIPE, stderr=STDOUT)
    if verbose or p.returncode != 0:
        print(p.stdout.decode('UTF-8'))
        
    if path.exists(outputdir + "/table.txt"):
        table = read_mumax3_table(outputdir + "/table.txt")
    else:
        table = None

    os.rename(scriptfile,outputdir+"/"+scriptfile)
        
    fields = read_mumax3_ovffiles(outputdir)
    
    return table, fields


def getScript(width,length,constant,spacing):

    resolution=2e-9
    constant=length/2*constant

    code="""
    resolution := """+str(resolution)+"""
    zResolution := 5e-09
    a := """+str(spacing)+"""

    gridSize := """+str(round(2*spacing/resolution))+""" //2*a/resolution
    gridDepth := 5 //25/zResolution

    SetCellsize(resolution,resolution,zResolution)
    SetGridsize(gridSize,gridSize,gridDepth)

    SetPBC(16, 16, 0)


    Msat = 700e3
    Aex = 13e-12
    alpha = 0.2
    Bmax := 0.1
    Bstep := Bmax / 200.0


    TableAdd(B_ext)
    TableAdd(m_full)


    islandWidth := """+str(width)+"""
    islandLength := """+str(length)+"""
    ellipseConstant:="""+str(constant)+"""


    hIsland := rect(islandLength-2*ellipseConstant, islandWidth)
    hIsland = hIsland.add(ellipse(ellipseConstant*2,islandWidth).transl(islandLength/2-ellipseConstant, 0, 0))
    hIsland = hIsland.add(ellipse(ellipseConstant*2,islandWidth).transl(-islandLength/2+ellipseConstant, 0, 0))
        

    vIsland := hIsland.rotz(pi / 2)
    vIslands := Universe().Inverse()

    hIslands := Universe().Inverse()
    hIslands = hIslands.add(hIsland.transl(-a/2, a, 0))
    hIslands = hIslands.add(hIsland.transl(-a/2, -a, 0))
    hIslands = hIslands.add(hIsland.transl(a/2, a, 0))
    hIslands = hIslands.add(hIsland.transl(a/2, -a, 0))
    hIslands = hIslands.add(hIsland.transl(-a/2, 0, 0))
    hIslands = hIslands.add(hIsland.transl(a/2, 0, 0))

    vIslands = vIslands.add(vIsland.transl(-a, a/2, 0))
    vIslands = vIslands.add(vIsland.transl(a, a/2, 0))
    vIslands = vIslands.add(vIsland.transl(-a, -a/2, 0))
    vIslands = vIslands.add(vIsland.transl(a, -a/2, 0))
    vIslands = vIslands.add(vIsland.transl(0, a/2, 0))
    vIslands = vIslands.add(vIsland.transl(0, -a/2, 0))

    DefRegion(1, hIslands)
    DefRegion(2, vIslands)
    TableAdd(m_full.Region(1))
    TableAdd(m_full.Region(2))
    setgeom(hIslands.add(vIslands))


    r2o2 := sqrt(2) / 2
    m = RandomMagSeed(4)
    B_ext = vector(Bmax*r2o2, Bmax*r2o2, 0)
    relax()
    minimize()

    for B := Bmax; B >= -Bmax; B -= Bstep {
        B_ext = vector(B*r2o2, B*r2o2, 0)
        minimize()
        tablesave()
    }
    for B := -Bmax; B <= Bmax; B += Bstep {
            B_ext = vector(B*r2o2, B*r2o2, 0)
            minimize()
            tablesave()
        }

    """
    
    return code

width=80e-9
#length=220e-9

pointyConstantMin=0
pointyConstantMax=1
pointyConstantStepCount=50

spacingMin=280e-9
spacingMax=500e-9
spacingStepCount=10

spacingVals=np.array([320])*1e-9

#2*spacing/resolution should have a lot of factors of 2

lengthVals=np.array([100,150,200,300,400])*1e-9
lengthMin=100e-9
lengthMax=500e-9
lengthStepCount=5


for constant in np.linspace(pointyConstantMin,pointyConstantMax,pointyConstantStepCount+1):
    for spacing in spacingVals:
        for length in lengthVals:
            if (spacing-length)/2<width/2:
                continue#islands have merged
            run_mumax3(getScript(width,length,constant,spacing), name=f"pointiness{constant}-spacing{spacing}-length{length}", verbose=False)
