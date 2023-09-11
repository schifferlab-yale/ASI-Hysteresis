import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import os
import glob
import shutil

mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['font.family'] = 'Helvetica'





def loadFile(name):
    data=pd.read_csv(name,sep="\t")
    data["H"]=np.round(data["B_extx (T)"]*10000/(np.sqrt(2)/2))
    data["m_full"]=(data["m_fullx (A/m)"]+data["m_fully (A/m)"])/np.sqrt(2)
    data["m"]=data["m_full"]/np.max(np.abs(data["m_full"]))

    try: 
        data["hIsland_my"]=data["m_full.region1y (A/m)"]/np.max(np.sqrt(data["m_full.region1y (A/m)"]**2+data["m_full.region1x (A/m)"]**2))
    except Exception:
        pass
    

    data=data.drop(columns=["# t (s)",\
    "mz ()","mx ()", "my ()",\
    "B_extx (T)", "B_exty (T)","B_extz (T)",\
    "m_fullz (A/m)","m_fullx (A/m)","m_fully (A/m)",\
    "m_full.region1x (A/m)","m_full.region1y (A/m)","m_full.region1z (A/m)",\
    "m_full.region2x (A/m)","m_full.region2y (A/m)","m_full.region2z (A/m)"])

    try:
        data.attrs["Hc"]=np.average(np.abs(data[np.sign(data["m"]).diff() != 0]["H"].iloc[1:]))
        data.attrs["Br"]=np.average(np.abs(data[np.sign(data["H"]).diff() != 0]["m"].iloc[1:]))
    except Exception:
        pass


    return data

def genMovie(path,name="out.mp4",pointiness=0,spacing=0,length=0):
    imagePaths=sorted(glob.glob(path+"/*.png"))[1:]
    images=[]
    data=loadFile(path+"/table.txt")

    for imPath in imagePaths:
        im = plt.imread(imPath)
        images.append(im)

    assert len(images)==len(data)



    #Make temporary directory
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, r'Temp Files')
    if os.path.exists(final_directory):
        shutil.rmtree(final_directory)
    os.makedirs(final_directory)


    for i in range(len(images)):

        #if i%200!=0: 
            #continue
            #print(i)

        fig,(imAx,pltAx)=plt.subplots(1,2)
        fig.set_size_inches(13, 7)
        fig.set_dpi(200)

        thisRow=data.loc[i,:]

        imAx.imshow(images[i],interpolation='nearest', aspect='auto')
        imAx.invert_yaxis()
        imAx.set_xlabel("x (nm)",fontsize=20)
        imAx.set_ylabel("y (nm)",fontsize=20)

        

        arrowScale=thisRow.H/max(data.H)

        imAx.arrow(0.5,0.5,0.3*arrowScale,0.3*arrowScale,width=0.1,head_length=0.1,lw=0,color="white",head_width=0.15,alpha=0.5,transform=imAx.transAxes)
        imAx.text(0.5,0.5,f"{round(thisRow.H)}Oe",horizontalalignment='center',
                verticalalignment='center',rotation=45,fontsize=20,fontweight=2,transform=imAx.transAxes)


        pltAx.plot(data["H"],data["m"],lw=3,color="k")
        pltAx.plot(-data["H"],-data["m"],lw=3,color="k")

        
        pltAx.scatter(thisRow["H"],thisRow["m"],s=200,zorder=10,color="red")
        pltAx.set_xlabel("H (Oe)", fontsize=20)
        pltAx.set_ylabel("$\\rm M$/$\\rm M_s$",fontsize=20)

        pltAx.tick_params(direction='in',top=1,right=1,width=1,length=4)

        plt.suptitle(f"p={round(pointiness*100)/100}\na={round(spacing*1e9)}nm\nl={round(length*1e9)}nm",fontsize=20)

        pltAx.tick_params(axis='x', labelsize=20)
        pltAx.tick_params(axis='y', labelsize=20)
        imAx.tick_params(axis='x', labelsize=20)
        imAx.tick_params(axis='y', labelsize=20)
        
        plt.tight_layout()
        plt.savefig(final_directory+f"/{i}.png")
        plt.close()

        os.system(f"ffmpeg -framerate 30  -i '{final_directory}/%d.png' \
  -c:v libx264 -pix_fmt yuv420p '{name}' -y >/dev/null 2>&1")





def decodeAttributes(string):
    dict={}
    desc=string.split(";")
    for attribute in desc:
        dict[attribute[0]] = float(attribute[1:])
    return dict["p"], dict["a"], dict["l"], dict["s"]


filePrefix="movieData/"
for name in glob.glob(filePrefix+"*.out")[0:10]:
    try:
        thisData=loadFile(name+"/table.txt")
    except Exception:
        print(f"Could not load {name}")
        continue
    desc=name[len(filePrefix):-len(".out")]
    print(name)
    print(desc)
    pointiness, spacing, length, runNum = decodeAttributes(desc)

    
    outFolder="movieOut"

    if not os.path.exists(os.path.join(os.getcwd(), outFolder)):
        os.makedirs(os.path.join(os.getcwd(), outFolder))

    outName=f"{outFolder}/p={round(pointiness*100)/100}; a={round(spacing*1e9)}nm; l={round(length*1e9)}nm"
    #outName=outName.replace(".","\.")
    outName+=".mp4"
    print(outName)

    genMovie(name,outName,pointiness,spacing,length)