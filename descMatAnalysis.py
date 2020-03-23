# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
import pandas as pd
from scipy.signal import lfilter
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyfit

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

def findYieldLimit(series):
    aux1 = series[:]


arq = pd.read_csv("a.csv", sep = ";")
arqB = pd.read_csv("b.csv", sep = ";")
arqC = pd.read_csv("c.csv", sep = ";")

arqB = arqB.apply(lambda x: np.float64(x))
arqC = arqC.apply(lambda x: np.float64(x))

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

# arq["Tension (GPa)"] /= 1e3
# arqB["Tension (GPa)"] /= 1e3
# arqC["Tension (GPa)"] /= 1e3

arq = arq.rename({"Force (N)":"Force (kN)", "Time (min)": "Time (s)", "Position (mm)":"Deformation (%)"}, axis = 1)
arqB = arqB.rename({"Position (mm)":"Deformation (%)"}, axis = 1)
arqC = arqC.rename({"Position (mm)":"Deformation (%)"}, axis = 1)
for name in ("Force (kN)","Strain (%)","Deformation (%)"):
    arq[name] = linFilter(arq[name])
    arqB[name] = linFilter(arqB[name])
    arqC[name] = linFilter(arqC[name])


# %%
start, end = 0.11, 0.35
force, strain = arq["Tension (GPa)"], arq["Strain (%)"]
boolArr = np.logical_and(strain >= start, strain <= end)
b, m = polyfit(strain[boolArr], force[boolArr], 1)
px, py = (0.2,0)
axis = np.linspace(0.2,0.65, 1000)
young_line = (axis-px)*m - py
young_mod, max_strain = intersectCurves(strain, force, axis, young_line)
plt.rcParams.update({'font.size': 22})
plt.figure(figsize = (21,15))
plt.scatter(young_mod, max_strain, c = "Blue", lw = 4, zorder = 5, label = "Módulo de Young")
plt.axhline(max_strain, ls = "--")
plt.plot(strain[boolArr], b + m*strain[boolArr], label = "Tendência", c = "green", lw = 3)
plt.plot(strain, force, label = "Amostra", alpha = 0.5, c = "gray")
plt.plot(axis, young_line)
plt.axvline(start, c = "red", ls = "--", label = "Intervalo Retilíneo")
plt.axvline(end,  c = "red", ls = "--")
plt.title("Curvas de Tensão da Amostra A")
plt.xlabel("Strain (%)")
plt.ylabel("Tension (GPa)")
# plt.xlim(0,1.09)
# plt.ylim(0,0.4)
# plt.annotate("%.2f" % max_strain, (0, max_strain + 0.5))
# plt.annotate("%.2f" % m, (0.8, 0.020))
plt.legend(loc = 4)
plt.savefig("Curvas de Tensão da Amostra A.png")
plt.show()
# print(m)


# %%
start, end = 0.03, 0.110
force, strain = arqB["Tension (GPa)"], arqB["Strain (%)"]
boolArr = np.logical_and(strain >= start, strain <= end)
b, m = polyfit(strain[boolArr], force[boolArr], 1)
px, py = (0.2,0)
axis = np.linspace(0.2,0.65, 1000)
young_line = (axis-px)*m - py
young_mod = (b + m*strain[boolArr]).mean()
plt.rcParams.update({'font.size': 22})
plt.figure(figsize = (21,15))
plt.plot(strain[boolArr], b + m*strain[boolArr], label = "Tendência", c = "green", ls = "-.", lw = 3)
plt.plot(strain, force, label = "Amostra", alpha = 0.5)
# plt.axhline(young_mod, lw = 3, label = ("Módulo de Young"))
plt.axvline(start, c = "red", ls = "--", label = "Intervalo Retilíneo")
plt.axvline(end,  c = "red", ls = "--")
plt.title("Curvas de Tensão da Amostra B")
plt.xlabel("Strain (%)")
plt.ylabel("Tension (GPa)")
# plt.xlim(0,1.09)
# plt.ylim(0,80)
plt.legend(loc = 4)
# plt.annotate("%.2f" % young_mod, (0, young_mod + 0.5))
# plt.annotate("%.2f" % m, (0.8, 30))
plt.savefig("Curvas de Tensão da Amostra B.png")
plt.show()
# print(m)


# %%
start, end = 0.001, 0.15
force, strain = arqC["Tension (GPa)"], arqC["Strain (%)"]
boolArr = np.logical_and(strain >= start, strain <= end)
b, m = polyfit(strain[boolArr], force[boolArr], 1)
px, py = (0.2,0)
axis = np.linspace(0.2,0.65, 1000)
young_line = (axis-px)*m - py
young_mod = (b + m*strain[boolArr]).mean()
plt.rcParams.update({'font.size': 22})
plt.figure(figsize = (21,15))
plt.plot(strain[boolArr], b + m*strain[boolArr], label = "Tendência", c = "green", ls = "-.", lw = 3)
plt.plot(strain, force, label = "Amostra", alpha = 0.5, c = "green")
plt.axhline(max(force), lw = 3, label = ("Módulo de Young"))
plt.axvline(start, c = "red", ls = "--", label = "Intervalo Retilíneo")
plt.axvline(end,  c = "red", ls = "--")
plt.title("Curvas de Tensão da Amostra C")
plt.xlabel("Strain (%)")
plt.ylabel("Tension (GPa)")
plt.xlim(0,1.09)
# plt.ylim(0,100)
plt.legend(loc = 4)
# plt.annotate("%.2f" % max(force), (0, max(force) + 0.5))
# plt.annotate("%.2f" % m, (0.8, 30))
plt.savefig("Curvas de Tensão da Amostra C.png")
plt.show()
# print(m)


# %%


