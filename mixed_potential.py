import numpy
import matplotlib.pyplot as plt 
import matplotlib.scale as mscale
from numpy import inf, log as ln, minimum, maximum
from numpy import log10 as log
from matplotlib.ticker import FixedLocator, NullFormatter
from numpy.core.function_base import geomspace
from scipy import interpolate
from scipy.optimize import fsolve

#given constants for simulation
R = 8.3144598
F = 96485
TafelC = 0.1
T = 298.15

#setting up files and graph
f = open('reactions.csv')
fig, (ax1) = plt.subplots(1)
ax1.set_xscale("log")
plt.xlabel("Current Density "+ "(A/cm^2)")
plt.ylabel("Potential (E(V) vs SHE)")
plt.title("Tafel Plot Slopes with Mixed Potentials")
mixedVoltageDataXOxi = []
mixedVoltageDataYOxi = []
mixedVoltageDataXRed = []
mixedVoltageDataYRed = []

#iterates through each reaction
for row in f:
    l = row
    name, standardPotential, Aconc, Acoeff, Bconc, Bcoeff, m, n, iLc, iLa, *_rest = l.split(',')
    standardPotential = float(standardPotential)
    Aconc = float(Aconc)
    Acoeff = float(Acoeff)
    Bconc = float(Bconc)
    Bcoeff = float(Bcoeff)
    m = int(m)
    n = int(n)
    iLc = float(iLc)
    iLa = float(iLa)

    #calculates Tafel plots for each reaction
    e = standardPotential + 2.3 * R * T / n / F * ln(Aconc**Acoeff / Bconc ** Bcoeff)
    ic = numpy.geomspace(Aconc, iLc, 10000)
    ia = numpy.geomspace(Bconc, iLa, 10000)
    yc = []
    ya = []
    maxi = e
    for x in ic:
        val = e + TafelC * log(x/Bconc) - 2.3*R*T/n/F * log(1-x/iLc)
        if (val == inf):
            yc.append(maxi)
        else:
            yc.append(val)
            maxi = val
    mini = e
    for x in ia:
        val = e - TafelC * log(x/Aconc) + 2.3*R*T/n/F * log(1-x/iLa)
        if (val == inf):
            ya.append(mini)
        else:
            ya.append(val)
            mini = val

    #adds to mixed potential reaction for oxidizing reaction
    if (len(mixedVoltageDataYOxi) == 0):
        for i in range(0, len(ic)):
            mixedVoltageDataXOxi.append(ic[i])
            mixedVoltageDataYOxi.append(yc[i])
    else:
        oxiFunction = interpolate.interp1d(mixedVoltageDataYOxi, mixedVoltageDataXOxi)
        newOFunction = interpolate.interp1d(yc, ic)

        rangeBelowDown = minimum(min(mixedVoltageDataYOxi), min(yc))
        rangeBelowUp = minimum(min(mixedVoltageDataYOxi), max(yc))
        rangeAboveDown = maximum(max(mixedVoltageDataYOxi), min(yc))
        rangeAboveUp = maximum(max(mixedVoltageDataYOxi), max(yc))

        #computes new oxidizing current densities and potentials for range of existing mixed voltage function
        for i in range(0, len(mixedVoltageDataYOxi)):
            try:
                mixedVoltageDataXOxi[i] = newOFunction(mixedVoltageDataYOxi[i]) + mixedVoltageDataXOxi[i]
            except ValueError:
                pass

        #computes new oxidizing current densities and potentials for values below existing range of existing mixed voltage function 
        try:
            for value in numpy.geomspace(rangeBelowUp, rangeBelowDown):
                mixedVoltageDataXOxi.insert(0, newOFunction(value))
                mixedVoltageDataYOxi.insert(0, value)
        except ValueError:
            pass

        #computes new oxidizing current densities and potentials for values above existing range of existing mixed voltage function
        try:
            for value in numpy.geomspace(rangeAboveDown, rangeAboveUp):
                mixedVoltageDataXOxi.append( newOFunction(value) + max(mixedVoltageDataXOxi))
                mixedVoltageDataYOxi.append(value)
        except ValueError:
            pass

    #adds to mixed potential reaction for reducing reaction 
    if (len(mixedVoltageDataYRed) == 0):
        for i in range(0, len(ia)):
            mixedVoltageDataXRed.append(ia[i])
            mixedVoltageDataYRed.append(ya[i])
    else:
        redFunction = interpolate.interp1d(mixedVoltageDataYRed, mixedVoltageDataXRed)
        newRFunction = interpolate.interp1d(ya, ia)

        rangeBelowDown = minimum(min(mixedVoltageDataYRed), min(ya))
        rangeBelowUp = minimum(min(mixedVoltageDataYRed), max(ya))
        rangeAboveDown = maximum(max(mixedVoltageDataYRed), min(ya))
        rangeAboveUp = maximum(max(mixedVoltageDataYRed), max(ya))

        #computes new reducing current densities and potentials for range of existing mixed voltage function
        for i in range(0, len(mixedVoltageDataYRed)):
            try:
                mixedVoltageDataXRed[i] = newRFunction(mixedVoltageDataYRed[i]) + mixedVoltageDataXRed[i]
            except ValueError:
                pass

        #computes new reducing current densities and potentials for values above existing range of existing mixed voltage function
        try:
            for value in numpy.geomspace(rangeAboveDown, rangeAboveUp):
                mixedVoltageDataXRed.insert(0, newRFunction(value))
                mixedVoltageDataYRed.insert(0, value)
        except ValueError:
            pass

        #computes new reducing current densities and potentials for values below existing range of existing mixed voltage function
        try:
            for value in numpy.geomspace(rangeBelowUp, rangeBelowDown):
                mixedVoltageDataXRed.append(newRFunction(value) + max(mixedVoltageDataXRed))
                mixedVoltageDataYRed.append(value)
        except ValueError:
            pass

    #graphs Tafel plots for individual reactions
    plt.plot(ic, yc, color = "black") 
    plt.plot(ia, ya, color = "black")
    #adds label to plot
    plt.annotate(name, xy=(ic[0], yc[0]), xytext=(ic[0], yc[0]+0.1), bbox = dict(boxstyle="round", fc="0.8"), arrowprops=dict(facecolor='gray', color = "gray", width = 0.01, headwidth = 0.01))

#graphs oxidizing and reduced reactions
plt.plot(mixedVoltageDataXOxi, mixedVoltageDataYOxi, color = "red", label = "mixed oxidizing reaction")
plt.plot(mixedVoltageDataXRed, mixedVoltageDataYRed, color = "blue", label = "mixed reduced reaction")
plt.legend()
plt.show()
