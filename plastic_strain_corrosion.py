import numpy
import matplotlib.pyplot as plt 
import matplotlib.scale as mscale
from numpy import inf, log as ln, minimum
from numpy import log10 as log
from matplotlib.ticker import FixedLocator, NullFormatter
from scipy.optimize import fsolve
import matplotx

#constants for Fe
P = 16666.6667
metal = "Fe"

#given constants for simulation
Vm = 9.99E-6
v = 0.45
a = 1.67E15
N0 = 1E12
R = 8.3144598
F = 96485
TafelC = 0.1
T = 298.15

#setting up files and graph
f = open('reactions_strain.csv')

fig, (ax1) = plt.subplots(1)
ax1.set_xscale("log")

plt.xlabel("Current Density") 
plt.ylabel("Voltage") 

reactants = {}

plt.xlabel("Current Density"+ "(A/cm^2)")
plt.ylabel("Potential (E(V) vs SHE)")
plt.title("Tafel Plot Slopes for H and " + metal)

#iterates through each reaction- calculating reaction for electrochemical only
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

    #calculating Tafel slopes for electrochemical reaction
    e = standardPotential + 2.3 * R * T / n / F * ln(Aconc**Acoeff / Bconc ** Bcoeff)
    reactants[name] = [e, Aconc, Acoeff, Bconc, Bcoeff, m, n, iLc, iLa]

    Ec = numpy.linspace(e, 1, 10000)
    Ea = numpy.linspace(e, -1, 10000)
    def oxi(Ec):
        return Aconc*10**((Ec-e)/TafelC)
    def red(Ea):
        return Bconc*10**((e-Ea)/TafelC)

    plt.plot(oxi(Ec), Ec, color = "black") 
    plt.plot(red(Ea), Ea, color = "black")
    plt.annotate(name, xy=(oxi(e), e+0.1), bbox = dict(boxstyle="round", fc="0.8"))
plt.show()

#factoring in stress corrosion
plt.xlabel("Plastic Strain (%)")
plt.ylabel("Potential (E(V) vs SHE)")
plt.title("Plastic Strain vs Corrosion Equilibrium Potential for " + metal)

if (reactants["H"][0] > reactants[metal][0]):
    for strain in numpy.geomspace(1, 1000):
        Ec = numpy.linspace(reactants[metal][0], 1, 10000)
        Ea = numpy.linspace(reactants["H"][0], -1, 10000)

        #incorporating continous elasto-plastic tension
        def oxi(Ec):
            return reactants["H"][1]*10**((Ec-reactants[metal][0]- P*Vm/(reactants[metal][6]*F)-T*R/(reactants[metal][6]*F)*ln(v*a/N0 * strain + 1))/TafelC)
        def red(Ea):
             return reactants["H"][1]*10**((reactants["H"][0]-Ea)/TafelC)*10**(strain*Vm/(-6*F*TafelC))
        def solve(E):
            return oxi(E)-red(E)

        plt.scatter(strain, fsolve(solve, 0.0001)[0])
else:
    for strain in numpy.geomspace(0.01, 0.4):
        Ec = numpy.linspace(reactants["H"][0], 1, 10000)
        Ea = numpy.linspace(reactants[metal][0], -1, 10000)

        #incorporating continous elasto-plastic tension
        def oxi(Ec):
            return reactants["H"][1]*10**((Ec-reactants["H"][0])/TafelC)
        def red(Ea):
             return reactants[metal][1]*10**((reactants[metal][0]-Ea- P*Vm/(reactants["H"][6]*F)-T*R/(reactants["H"][6]*F)*ln(v*a/N0 * strain + 1))/TafelC)*10**(strain*Vm/(-6*F*TafelC))
        def solve(E):
            return oxi(E)-red(E)

        plt.scatter(strain, fsolve(solve, 0.0001)[0], color = "blue")
plt.show()
