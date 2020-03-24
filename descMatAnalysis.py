# De autoria de Carlos Eduardo Dip, João Pedro Gianfaldoni de Andrade e Lucas Keichi Fukada.

print("Importing...")
import pandas as pd
from scipy.signal import lfilter
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyfit
import os
print("Import Successful!")
print("Reading Files...")
plt.rcParams.update({'font.size': 22})
    
arq = pd.read_csv(r"C://Users//cadud//Repos//MaterialAnalysis//data//a.csv", sep = ";")
arqB = pd.read_csv(r"C://Users//cadud//Repos//MaterialAnalysis//data//b.csv", sep = ";")
arqC = pd.read_csv(r"C://Users//cadud//Repos//MaterialAnalysis//data//c.csv", sep = ";")
# print(arq)
arqB = arqB.apply(lambda x: np.float64(x))
arqC = arqC.apply(lambda x: np.float64(x))
print("Files read successfully!")
print("Correcting units and sorting data...")
# arq["Force (N)"] /= 1000
arq["Time (min)"] *= 60
# arq["Position (mm)"] /= 50
# arqB["Position (mm)"] /= 50
# arqC["Position (mm)"] /= 50
arqB["Force (kN)"] *= 1000
arqC["Force (kN)"] *= 1000

arq["Tension (GPa)"] = arq["Force (N)"]/(np.pi*(8.5/2)**2) #arq["Force (N)"]/(np.pi*(8.5/2000)**2)
arqB["Tension (GPa)"] = arqB["Force (kN)"]/(np.pi*(8.55/2)**2)
arqC["Tension (GPa)"] = arqC["Force (kN)"]/(np.pi*(8.5/2)**2) 

arq["Tension (GPa)"] /= 10
arqB["Tension (GPa)"] /= 10
arqC["Tension (GPa)"] /= 10

arq = arq.rename({"Force (N)":"Force (kN)", "Time (min)": "Time (s)"}, axis = 1)
# arqB = arqB.rename({"Position (mm)":"Deformation (%)"}, axis = 1)
# arqC = arqC.rename({"Position (mm)":"Deformation (%)"}, axis = 1)
print("Data organized successfully!")
def linFilter(y):
    n = 55  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    yy = lfilter(b,a,y)
    return yy

def intersectCurves(c1x, c1y, c2x, c2y):
    force = c1y
    strain = c1x
    axis = c2x
    young_line = c2y
    lista = []
    for y1, x1 in zip(force, strain):
        for x2, y2 in zip(axis, young_line):
            if np.allclose(np.array([x1, y1]), np.array([x2,y2]), 0.01):
                lista.append([x2,y2])
    lista = np.array(lista).mean(axis = 0)
    return (lista[0], lista[1])

def findLocalMaximum(strain, start, end, force):
    return max(force[np.logical_and(strain>start,strain<end)])

def Plot(start = 0.11, end = 0.35, force = arq["Tension (GPa)"], strain = arq["Strain (%)"], draw_projection = False, writings = {"title":"a", "xlabel":"b", "ylabel":"c", "filename":"unnamed.png"}, xbounds = (0,1), ybounds = (0,30), drawMaxStrain = False, limitsForLocalMax = None):
    print("Started drawing curve %s" % writings["title"])
    
    if not writings["title"]:
        writings["title"] = "a"
    if not writings["xlabel"]:
        writings["xlabel"] = "b"
    if not writings["ylabel"]:
        writings["ylabel"] = "c"
    if not writings["filename"]:
        writings["filename"] = "unnamed.png"
    
    ## Fazer o gráfico.
    plt.figure(figsize = (21,15))

    if ((start is not None) and (end is not None)):
        print("Drawing trendline...")
        # Determina a região a ser considerada reta.
        boolArr = np.logical_and(strain >= start, strain <= end)
        
        # Faz um ajuste de retas por MMQ para a região acima.
        b, youngMod = polyfit(strain[boolArr], force[boolArr], 1)
            
        # Marca o intervalo para análise
        plt.axvline(start, c = "red", ls = "--", label = "Intervalo Retilíneo")
        plt.axvline(end,  c = "red", ls = "--")

        # Desenha a linha de tendência.
        plt.plot(strain[boolArr], b + youngMod*strain[boolArr], label = "Tendência", c = "green", lw = 3)
        
        # Escreve o módulo de Young no gráfico.
        plt.scatter([0],[0], label = "Módulo de Young: %.2f (GPa)" % youngMod, alpha  = 0)
        print("Trendline drawn!")


    # Determina a região a ser olhada para interseção
    boolArr2 = np.logical_and(strain >= end, strain <= 1)

    # Determina a região a ser olhada para determinar o limite de resistência.

    if limitsForLocalMax is not None:
        _limsmall, _limlarge = limitsForLocalMax
        _max = findLocalMaximum(strain, _limsmall, _limlarge, force)
        plt.axhline(_max, label = "Limite de Escoamento: %.2f (GPa)" % _max, alpha = 0.7, c = "Green", ls = "--")


    if draw_projection:
        print("Searching for intersection in curve %s..." % writings["title"])
        px, py = (0.2,0)
        axis = np.linspace(0.2,0.65, len(strain[boolArr2]))
        young_line = (axis-px)*youngMod - py
        elasticLimit, intersectedStrain = intersectCurves(strain[boolArr2], force[boolArr2], axis, young_line)
        print("Intersection Found!")
        plt.scatter(elasticLimit, intersectedStrain, c = "Blue", lw = 4, zorder = 5, label = "Limite de Escoamento: %.2f (GPa)" % intersectedStrain)
        plt.plot(axis, young_line)
        plt.annotate("%.2f" % intersectedStrain, (elasticLimit + 0.5, intersectedStrain + 0.5))


    
    print("Drawing your graph now!")
    if drawMaxStrain:
        # Desenha uma linha horizontal no limite de resistência.
        plt.axhline(max(force), label = "Limite de Resistência: %.2f (GPa)"%max(force), alpha = 0.7, c = "Green", ls = "--")

    # Desenha a curva de deformação por tensão.
    plt.plot(strain, force, label = "Amostra", alpha = 0.5, c = "black")
    plt.title(writings["title"])
    plt.xlabel(writings["xlabel"])
    plt.ylabel(writings["ylabel"])
    if xbounds is not None:
        plt.xlim(xbounds[0], xbounds[1])
    if ybounds is not None:    
        plt.ylim(ybounds[0],ybounds[1])
    plt.legend(loc = 4)
    plt.savefig(writings["filename"])
    plt.show()

print("Filtering data to remove bumps...")
for name in ("Force (kN)","Strain (%)","Position (mm)"):
    arq[name] = linFilter(arq[name])
    arqB[name] = linFilter(arqB[name])
    arqC[name] = linFilter(arqC[name])
print("Data filtered succesfully!")

### ----------------------------- ###
### Começo das curvas individuais ###
doA = True
doB = True
doC = True

if doA:
    writings = {"title":"Amostra A",
                "xlabel":"Deformação relativa (%)", 
                "ylabel":"Tensão (GPa)", 
                "filename":"Str_Tens(A).png"}
    Plot(0.11, 0.35, arq["Tension (GPa)"], arq["Strain (%)"], draw_projection=True, writings = writings, ybounds = (0,31))

    ### ----

    writings = {"title":"Amostra A",
                "xlabel":"Deformação absoluta (mm)", 
                "ylabel":"Tensão (GPa)", 
                "filename":"Disp_Tens(A).png"}
    Plot(None, None, arq["Tension (GPa)"], arq["Position (mm)"], draw_projection=False, writings = writings, xbounds = (0,9), ybounds = (0,35), drawMaxStrain=True)

### ----

if doB:
    writings = {"title":"Amostra B",
                "xlabel":"Deformação relativa (%)", 
                "ylabel":"Tensão (GPa)", 
                "filename":"Str_Tens(B).png"}
    Plot(0.03, 0.11, arqB["Tension (GPa)"], arqB["Strain (%)"], draw_projection=False, writings = writings, ybounds = (0,30), limitsForLocalMax = (0.01, 0.4))

    ### ----

    writings = {"title":"Amostra B",
                "xlabel":"Deformação absoluta (mm)", 
                "ylabel":"Tensão (GPa)", 
                "filename":"Disp_Tens(B).png"}
    Plot(None, None, arqB["Tension (GPa)"], arqB["Position (mm)"], draw_projection=False, writings = writings, xbounds = (0,25), ybounds = (0,41), drawMaxStrain = True)

### ----

if doC:
    writings = {"title":"Amostra C",
                "xlabel":"Deformação relativa (%)", 
                "ylabel":"Tensão (GPa)", 
                "filename":"Str_Tens(C).png"}
    Plot(0.001, 0.15, arqC["Tension (GPa)"], arqC["Strain (%)"], draw_projection=False, writings = writings, ybounds = (0,41), limitsForLocalMax = (0.15, 0.25))

    ### ----

    writings = {"title":"Amostra C",
                "xlabel":"Deformação absoluta (mm)", 
                "ylabel":"Tensão (GPa)", 
                "filename":"Disp_Tens(C).png"}
    Plot(None, None, arqC["Tension (GPa)"], arqC["Position (mm)"], draw_projection=False, writings = writings, xbounds = (0,19), ybounds = (0,65), drawMaxStrain = True)