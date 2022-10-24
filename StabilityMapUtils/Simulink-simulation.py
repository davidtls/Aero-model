"""
Created on Friday Oct 21 21:35:26 2022

@author: David Planas Andrés
david.planas-andres@isae-supaero.fr



File for performing dynamic simulation with Simulink and Distributed Propulsion Wing-Propeller interaction.


The idea is:

 1) Defining here the parameters for the simulink simulation: testing time, time step...
 2) Defining here the initial flight situation
 3) Calling to the objects geometry and propwing and passing them to simulink
 4) Defining the initial values of : V, alpha, beta, gamma, theta, omega, p, q, r, phi
   delta_e, delta_a, delta_r and passing them to simulink

   Ideally the system should not be excited in the beggining ?


 5) Running the simulink directly from here with the help of the API. This is done with an API
    called matlab engine.


  The complicate thing is that one modulus of simulink must call Python:

  This does not look like complicate, just use a matlub function in simulink, that is an editor-like
  script of matlab where you can call python using the API:

     output1, output2 ... = py.Function (input1, input2, ...)




NOTES:
    Ejecutar en la línea de comandos


cd "C:\Program Files\MATLAB\R2020b\extern\engines\python"

python setup.py install

NO PUEDO PORQUE NO ME DEJA SALIR DE LA CARPETA C


import matlab.engine

eng = matlab.engine.start_matlab()     arranca matlab para poder hacerle llamadas

a partir de ahora pones eng. y el comando que sea de Matlab y lo tienes

"""







import numpy as np
import math
import scipy.linalg
import scipy.io  # input/output with matlab
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.optimize import least_squares
import sys
import PattersonAugmented as PA
import control
import pylab
import ReadFileUtils
from StabilityMapUtils import equation as e
from StabilityMapUtils import AeroForces
import time
import pickle
import datetime






"""
Defining the geometry g (see if some other important atributes are defined later) and H0, V0, gamma0, beta0, omega0
"""

from AircraftClass import ATRgeometry

Neng = 12
inop_eng = 0
FlapDefl = 0 * np.pi / 180  # in degree standard flap deflection. Deflections allowed : 0 15 and 30 degree. Keep in mind they can be deflected for V<=71

g = ATRgeometry.data(1.0, Neng, inop_eng, FlapDefl, TipClearance=True, dprop=0.1,
                     dfus=0.1)  # arg = Vtsize + options(Neng, inop_eng, vertical tail parameters)


# --- Dictionnary for type of aircraft studied. aircraft: ATR72, version : 'original', 'DEPoriginal', 'DEPnofin'
g.hangar = {'aircraft': 'ATR72', 'version': 'original'}

# --- Test case and steady parameters
Vsr = 50.9  # m/s at 21.5T with 15°Fl, or 59.2ms at 0°Fl
H_base = 0  # in m the altitude
V_base = 70          #1.3 * Vsr
beta_base = 0 / 180 * math.pi
gamma = 0
R = 000  # in meters the turn radius
phimax = 5  # in degree the max bank angle authorized
alphamax = 25  # in degree, stall bound for trimming
deltaRmax = 30  # in degree
ThrottleMax = 1  # max thrust level
ThrottleMin = 1e-9  # min throttle, don't accept 0 thrust
g.VelFlap = 0  # modificada para que no salte   71  # in m/s the maximum velocity at which flap are deployed





if g.IsPropWing:
    g.alpha_max = 11.7 / 180 * np.pi + g.alpha_i - g.alpha_0
    g.alpha_max_fl = 11.7 / 180 * np.pi + g.alpha_i - g.alpha_0



Velocities = (70, 90, 110, 130, 150)                          #Two first speeds for H = 0 m, last 3 for H=5000m
rho_vec = (1.225, 1.225, 0.736116,  0.736116,  0.736116)      #rho= 0.736116 kg/m^3     H=5000 m   a=320.529 m/s
Mach = [0.2058, 0.2647, 0.3431, 0.4055, 0.4679]







path = 'ATR72_SI_MTOW_FinLess_STAB/'
filenameNoFin = [path+'ATR72_FinLess_mach1.stab', path+'ATR72_FinLess_mach2.stab', path+'ATR72_FinLess_mach3.stab',
                        path +'ATR72_FinLess_mach4.stab', path+'ATR72_FinLess_mach5.stab']
Matrix = ReadFileUtils.ReadStabCoef(filenameNoFin)



CoefMatrix = g.NicolosiCoef(Matrix[:, 1:], Mach)



atmospher = g.GetAtmo(H_base)
a_sound = atmospher[0]
rho_base = atmospher[1]
M_base = V_base/a_sound





Coef_base = AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)
g.Matrix_no_tail_terms = AeroForces.CoefInterpol(M_base, Matrix[:, 1:], Mach)






PropPath = "./ATR72_SI_MTOW_Control_FinLess_FEM/"
PropFilenames = {'fem': [PropPath+"ATR72_FinLess_mach1",
                         PropPath+"ATR72_FinLess_mach2",
                         PropPath+"ATR72_FinLess_mach3",
                         PropPath+"ATR72_FinLess_mach4",
                         PropPath+"ATR72_FinLess_mach5"],
                 'AirfoilPolar': PropPath+"naca3318Pol.txt",
                 'FlapPolar': PropPath+"naca3318fl+10.txt",
                 'AileronPolar': PropPath+"naca3318fl+10.txt"}
PW = PA.PropWing(g, PropFilenames)
PW.DeltaCL_a_0 = 1

