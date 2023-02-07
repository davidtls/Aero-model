# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:05:08 2017

Main file for stability mapping

@authors: e.nguyen-van ,david.planas-andres
"""

# import all modules


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
import Eigenvalues_manually
import derivatives_calc
import Apricott
import Python_to_Matlab
import CL_figure


import Forces_test
from StabilityMapUtils import Longitudinal
from StabilityMapUtils import Equilibrated_polar as EQ


"""
Things to check:
         Aircraft model
         Propeller Wing Interaction
         inop engines
         Speed
         Altitude
         Complete Optimization/ Long_equilibrium
         Forces comparison deactivated
         g.hangar DEPoriginal (different thrust) or original (same thrust)
         flaps
         Velocity of flaps deployment
"""



aircraft = {'model': 'X-57'}   #OPTIONS: ATR, DECOL, X-57







if aircraft['model'] == 'DECOL':
    from AircraftClass import DECOLgeometry


    Neng = 8
    inop_eng = 0
    FlapDefl = 0 * np.pi / 180  # in degree standard flap deflection. Deflections allowed : 0 and 15  degree.
                                # For flap deflection V_base < VelFlap


    g = DECOLgeometry.data(1, Neng, inop_eng, FlapDefl , r=0.113 / 2, rf=0.1865 / 2, zw=0.045,TipClearance=True, dprop=0.1)  # arg = Vtsize + options(Neng, inop_eng, vertical tail parameters)

    # --- Test case and steady parameters
    H_base = 0  # in m the altitude
    V_base = 15  # 16.38
    beta_base = 0 / 180 * math.pi
    gamma = 0 / 180 * np.pi  # math.atan(0/87.4)#/180*math.pi # 3% slope gradient # 6.88m/s vertical
    R = 0  # in meters the turn radius
    phimax = 10  # in degree the max bank angle authorized
    alphamax = 25  # in degree to adjust if it is below 71m/s
    deltaRmax = 30  # in degree
    ThrottleMax = 1  # max thrust level
    ThrottleMin = 0.0001  # -0.34





    g.hangar = {'aircraft':'DECOL', 'version':'original'}  # Original for twin engine model (all engines same power)

    g.P_var = 8 * 14.4 * 4  # I*V*N_eng/2
    g.VelFlap = 12.5  # in m/s the maximum velocity at which flap are deployed
    g.alpha_max = 10 / 180 * np.pi
    g.alpha_max_fl = 10 / 180 * np.pi




    # FLight measured Cd0:

    # --- additional parameters (default edited during execution) ---
    g.set_nofin(False)  # =True means : no rudder used
    g.Pkeyword="Selig"



    # ---- Optim parameter ------
    MaxIter = 100
    tolerance = 1e-3
    method = 'SLSQP'





elif aircraft['model'] == 'ATR':
    from AircraftClass import ATRgeometry


    Neng = 12
    inop_eng = 0
    FlapDefl = 0 * np.pi / 180  # in degree standard flap deflection. Deflections allowed : 0 15 and 30 degree. Keep in mind they can be deflected for V<=71

    g = ATRgeometry.data(1.0, Neng, inop_eng, FlapDefl, TipClearance=True, dprop=0.1,
                         dfus=0.1)  # arg = Vtsize + options(Neng, inop_eng, vertical tail parameters)

    # --- Test case and steady parameters
    Vsr = 50.9  # m/s at 21.5T with 15°Fl, or 59.2ms at 0°Fl






    H_base = 0  # in m the altitude
    V_base = 70          #1.3 * Vsr
    beta_base = 0 / 180 * math.pi
    gamma = 0
    # Climb condition  np.arctan(3 / 100)  # (3/100)/180*np.pi##math.atan(0/87.4)#/180*math.pi # 3% slope gradient # Best climb rate: 6.88m/s vertical @ 87.5m/s = 4.5°gamma, see http://www.atraircraft.com/products_app/media/pdf/Fiche_72-600_Juin-2014.pdf
    R = 000  # in meters the turn radius
    phimax = 5  # in degree the max bank angle authorized
    alphamax = 25  # in degree, stall bound for trimming
    deltaRmax = 30  # in degree
    ThrottleMax = 1  # max thrust level
    ThrottleMin = 1e-9  # min throttle, don't accept 0 thrust
    ThrottleMinExt = 1e-9  # -0.34


    # --- dictionnary for type of aircraft studied. aircraft: ATR72, version : 'original', 'DEPoriginal', 'DEPnofin'
    g.hangar = {'aircraft': 'ATR72', 'version': 'original'}

    g.VelFlap = 0  # modificada para que no salte   71  # in m/s the maximum velocity at which flap are deployed
    g.CL0_fl = g.CL0_fl  # / 2   # for take off in no-interaction divide by two








    if g.IsPropWing:
        # ensures alpha fuselage as reference for stall
        g.alpha_max = 11.7 / 180 * np.pi + g.alpha_i - g.alpha_0
        g.alpha_max_fl = 11.7 / 180 * np.pi + g.alpha_i - g.alpha_0

        #Explanation
        """
        The idea is that alpha_max =  (alpha-stall-airfoil) + (alpha_0_airfoil) + offset(wing-airfoil).
        Effectively 11.7 + 4 = 15.7 is the angle of stall of the airfoil. We are also adding alpha_0 (alpha_0 is negative and we are
        substracting, so adding something at the end)  as we are comparing alpha_max with  (aoa) + alpha_i + alpha_0
        to see if there is stall or not. So at the end, all this is meaning that if the angle between the airspeed and the
        fusselage is 11.7, the airfoil has an incidence of 15.7 and is in stall. 
        """



    else:
        g.alpha_max = 11.7 / 180 * np.pi          # In Aeroforces the stall is not considered when there is no interaction
        g.alpha_max_fl = 11.7 / 180 * np.pi       # so not really a lot of sense in this section

    # --- additional parameters (default edited during execution) ---
    g.set_nofin(False)  # =    True means : no rudder used    False:rudder used




    """ Algorithm set up """
    # ---- Optim parameter ------
    MaxIter = 100
    tolerance = 1e-5
    method = 'SLSQP'  # 'trust-interior' or 'SLSQP'








elif aircraft['model'] == 'X-57':
    from AircraftClass import X57geometry
    from StabilityMapUtils import Longitudinal

    HLP = True
    if HLP:
        Neng = 12
        FlapDefl = 30 * np.pi / 180  # in degree standard flap deflection. Deflections allowed : 0 10 and 30 degree
        V_base = 40.145
        H_base = 762   # in m the altitude : 2434 m (cruise) / 0 m take-off
    else:
        Neng = 2
        FlapDefl = 0 * np.pi / 180
        V_base = 77.1677         # 77.1677 (cruise)   or stall speed (Vsr=29.84)
        H_base = 762  # in m the altitude : 2434 m (cruise) / 0 m take-off

    inop_eng = 0


    g = X57geometry.data(1.0, Neng, inop_eng, FlapDefl, HLP)

    # --- Test case and steady parameters
    Vsr = 29.84  # m/s with 30° flaps and high lift propellers




    beta_base = 0 / 180 * math.pi
    gamma = 0                                                                                                           # previous condition  np.arctan(3 / 100)  # (3/100)/180*np.pi##math.atan(0/87.4)#/180*math.pi # 3% slope gradient # Best climb rate: 6.88m/s vertical @ 87.5m/s = 4.5°gamma, see http://www.atraircraft.com/products_app/media/pdf/Fiche_72-600_Juin-2014.pdf
    R = 000  # in meters the turn radius
    phimax = 5  # in degree the max bank angle authorized
    alphamax = 25  # in degree, stall bound for trimming
    deltaRmax = 30  # in degree
    ThrottleMax = 1  # max thrust level
    ThrottleMin = 1e-9  # min throttle, don't accept 0 thrust
    ThrottleMinExt = 1e-9  # -0.34


    # --- dictionnary for type of aircraft studied. aircraft: ATR72, version : 'original', 'DEPoriginal', 'DEPnofin'
    g.hangar = {'aircraft': 'X-57', 'version': 'original'}

    g.VelFlap = 0  # 0, we will set it.
    g.CL0_fl = g.CL0_fl  # / 2   # for take off in no-interaction divide by two

    if HLP:
        g.IsPropWing = True
        g.IsPropWingDrag = True

        # ensures alpha fuselage as reference for stall. Stall angle from GNEW5BP93B airfoil documentation.
        g.alpha_max = 18 / 180 * np.pi - g.alpha_0    # See explanation for ATR and keep in mind here alpha_i is 0 (no wing-fuselage offset)
        g.alpha_max_fl = 12 / 180 * np.pi - g.alpha_0

    else:
        g.IsPropWing = False
        g.IsPropWingDrag = False







    # --- additional parameters (default edited during execution) ---
    g.set_nofin(False)  # =    True means : no rudder used    False:rudder used





    """ Algorithm set up """
    # ---- Optim parameter ------
    MaxIter = 100
    tolerance = 1e-5
    method = 'SLSQP'  # 'trust-interior' or 'SLSQP'





































mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif'] = ['Latin Modern Roman']  # for legend / FreeSerif
mpl.rcParams["mathtext.fontset"] = "stix"  # for math in legend
mpl.rcParams["font.size"] = 16 






print("Engine position :")
strout = ""
for i in range(len(g.yp)):
    strout = strout+str(i+1)+" : "+"{:.3}".format(g.yp[i])+"m, "
print(strout)





#Options
g.Complete_Optimization = False
g.Long_equilibrium = True





#Prop-wing interaction settings
g.DisplayPatterInfo = False





# --- Study jacobian
gojac = True
storeJac = False
OutputMatlab = True
goFinVariation = False
FinRatioVec = np.linspace(0.1, 1.0, 10)
CstSpan = False
CstA = True



















    






# --- hard coded velocity, corresponding air density and Sound vel ---
# The velocities corresponds to key points of ATR flights and the ones at which
# VSPaero computes the stability matrices

if aircraft['model'] == 'ATR':

     Velocities = (70, 90, 110, 130, 150)                                          #Two first speeds for H = 0 m, last 3 for H=5000m
     rho_vec = (1.225, 1.225, 0.736116,  0.736116,  0.736116)                      #rho= 0.736116 kg/m^3     H=5000 m   a=320.529 m/s
     # a_sound = (340, 340, 320.529 ,320.529, 320.529)
     Mach = [0.2058, 0.2647, 0.3431, 0.4055, 0.4679]                               # Mach= Vel/a


elif aircraft['model'] =='DECOL':

     Velocities = (10, 15, 20, 25, 30, 35)
     rho_vec = (1.225, 1.225, 1.225, 1.225, 1.225, 1.225)
     Mach = np.ones((len(Velocities), 1))*0.0001

elif aircraft['model'] =='X-57':

     Velocities = (30, 35, 40, 75, 80)
     rho_vec = (1.225, 1.225, 0.962870, 0.962870, 0.962870)
     Mach = [0.0882, 0.1029, 0.1175, 0.2267, 0.2418]











""" Algorithm start """

#--- List all .stab file from vsp aero and read the coeff ---- 

Matrix = ReadFileUtils.ReadStabCoef(g.filenameNoFin)
        

print("Adjusting Kf, Kh and VT size")
print("New Kf and Kh")
print(g.AdjustVT())

CoefMatrix = g.NicolosiCoef(Matrix[:, 1:], Mach)
        

# --- Interpol coefficients for test case ---
# Find sound velocity and air density
atmospher = g.GetAtmo(H_base)
a_sound = atmospher[0]
rho_base = atmospher[1]
M_base = V_base/a_sound



Coef_base = AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)
g.Matrix_no_tail_terms = AeroForces.CoefInterpol(M_base, Matrix[:, 1:], Mach)







# Define here the PropWing interaction

PW = PA.PropWing(g, g.PropFilenames)
PW.DeltaCL_a_0 = 1



if aircraft['model']=='DECOL':

            #PW.AoAZero[:,-1] = PW.AoAZero[:,-1] + 3.2/180*np.pi #correction for angle of incidence of wing
            PW.AoAZero[:, 0] = PW.AoAZero[:, 0]*10**(-3)
            PW.CLslope[:, 0] = PW.CLslope[:, 0]*10**(-3)
            PW.AoAZero[:, 1] = PW.AoAZero[:, 1]*10**(-6)
            PW.CLslope[:, 1] = PW.CLslope[:, 1]*10**(-6)         #TO CORRECT UNITS AS DECOL.vsp3 is in mm
            PW.AoAZero[:, 2] = PW.AoAZero[:, 2]*10**(-3)
            PW.CLslope[:, 2] = PW.CLslope[:, 2]*10**(-3)
            PW.DeltaCL_a_0 = 1  # CL_alpha correction factor













#Forces_comparison = Forces_test.Constraints_DEP(Coef_base, atmospher, g, PW)
#Forces_comparison = Forces_test.Constraints_DEP_body(Coef_base, atmospher, g, PW)
Forces_comparison = Forces_test.Long_equilibrium2(Coef_base, atmospher, g, PW)


# k, x, fixtest, diccons = EQ.Equilpolar(Matrix, CoefMatrix, Mach, g, PW)









if g.Complete_Optimization == True:


        #Initialise test and guess vectors
        if g.nofin == False:
            # x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
            x0 = np.array([5*math.pi/180, 0, 0, 0, 0.00, 0.0, 0.0, 0.0, 0])
            bnds = ((-5*math.pi/180, alphamax*math.pi/180), (-0.2, 0.2), (-0.2, 0.2), (-0.2, 0.2), (-phimax/180*math.pi, phimax/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-23/180*math.pi, 13/180*math.pi), (-deltaRmax/180*math.pi, deltaRmax/180*math.pi))
        else:
            # x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_i]
            x0 = np.array([7.5*math.pi/180, 0, 0, 0, 0.00, 0.05, 0.00, 0.00])
            bnds = ((-5*math.pi/180, alphamax*math.pi/180), (-0.2, 0.2), (-0.2, 0.2), (-0.2, 0.2), (-phimax/180*math.pi, phimax/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-23/180*math.pi, 13/180*math.pi))

        eng_vec = np.array([0.1] * g.N_eng)



        # complete the vectors with engines:
        x0 = np.append(x0, eng_vec)

        ## General formulation:
        bnds_eng = ((ThrottleMin, ThrottleMax), (ThrottleMin, ThrottleMax))
        for i in range(int(g.N_eng / 2)):
            bnds = bnds + bnds_eng

        # --- imposed conditions ---
        # fix = [V, beta, gamma, omega, H]
        if R == 0:
            omega = 0
        else:
            omega = V_base * math.cos(gamma) / R

        fixtest = np.array([V_base, beta_base, gamma, omega])

        # put everything in tuples for passing to functions
        diccons = (np.copy(fixtest), np.copy(Coef_base), atmospher, g, PW)  # fix, CoefMatrix,Velocities, rho, g
        dicfobj = (np.copy(fixtest), rho_base, g)



        # --- minimization algorithm ---
        ## SLSQP
        if method == 'SLSQP':
            k = minimize(e.fobjectivePropWingInterac, np.copy(x0), args=dicfobj, bounds=bnds,
                         constraints={'type': 'eq', 'fun': e.Constraints_DEP, 'args': diccons},
                         options={'maxiter': MaxIter, 'disp': True}, tol=tolerance)
        elif method == 'trust-interior':
            ## Trust interior
            k = minimize(e.fobjectivePropWingInterac, np.copy(x0), args=dicfobj, method='trust-constr', bounds=bnds,
                         constraints={'type': 'eq', 'fun': e.Constraints_DEP, 'args': diccons},
                         options={'maxiter': MaxIter, 'disp': True}, tol=tolerance)
        else:
            sys.exit('Non valid optimization method')

        x = k.x







if g.Long_equilibrium:

    """
     Function for longitudinal analysis.
     * 6 variables can be played with : alpha, gamma, V, de, dx, theta
     * There are a total of 4 equations, two for forces in axes x and z, one for moments, and one relating gamma,
       alpha and theta.

      Function has two ways of working:
     * If two variables are defined among these 6, then the function will give away the value of the other four.
     * If zero, or one variable are defined, then the function will give a set of 6 variables that optimize
       a given function that must be defined.
     * If more than two variables are defined problem may be over-constrained and no solution will be given.


     For defining a variable just be careful to give it a value inside the possible limits.
     For leaving a variable undefined just define it as a string.

    Examples:
    * By just defining gamma = 0 the function will find the minimum speed (stall speed) at the altitude given (while the
      objective function is V_min) (OPTIMIZATION PERFORMED).
    * By defining V and dx    the function will find the set of (alpha de gamma theta) to accomplish equilibrium.
    * By defining V and gamma the function will find the set of (alpha de dx    theta) to accomplish equilibrium.
    """

    #Long_variables;   (to define a combination of 2 or less)
    alpha_long = "TS"
    gamma_long = 0
    V_long = 40.145
    de_long = "TS"
    dx_long = "TS"
    theta_long = "TS"

    vars = [alpha_long, gamma_long, V_long, de_long, dx_long, theta_long]

    k, x, fixtest, diccons = Longitudinal.Long_Equilibrium(Coef_base, atmospher, g, PW, vars)













# print results
print(k)
def printx(x, fix, atmo, g, PW):
    V = fix[0]
    alpha = x[0]/math.pi*180
    beta = fix[1]/math.pi*180
    pqr = x[1:4]/math.pi*180
    phi = x[4]/math.pi*180
    theta = x[5]/math.pi*180
    da = x[6]/math.pi*180
    de = x[7]/math.pi*180

    print("\nState vector value:")
    print("V= {0:0.2f}m/s, alpha = {1:0.2f}\xb0, beta={2:0.2f}\xb0, phi={3:0.2f}\xb0, theta={4:0.2f}\xb0".format(V, alpha, beta, phi, theta))
    print("p={0:0.4f}\xb0/s q={1:0.4f}\xb0/s r={2:0.4f}\xb0/s".format(*pqr))
    print("da={0:0.2f}\xb0, de= {1:0.2f}\xb0".format(da, de))

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * fix[1] + g.wingsweep) - x[3] * g.yp

    if g.IsPropWing:
        if V <= g.VelFlap or g.FlapDefl != 0:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect, atmo)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], g.FlapDefl, g, False, fix[1], x[1], V, x[3])
        else:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect, atmo)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], 0, g, False, fix[1], x[1], V, x[3])

    if g.nofin == False:
        print("dr = {0:0.2f}\xb0".format(x[8]/math.pi*180))




printx(x, fixtest, atmospher, g, PW)















# check if constraints are validated
constraints_calc = e.Constraints_DEP(x, *diccons)
print("\nConstraints")
print(constraints_calc)

















# ---- Compute jacobian -----
if gojac == True:
    jac = e.Jac_DEP(x, *diccons, 0.01)


    #jacCol=[ V, beta, gamma, alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
    #jacLign = [V, beta, alpha, p,q,r,phi,theta)                                                                        element 2,3   (count starting on 1)     is      d beta /  d gamma

    #Add a lign to jac : gamma_dot = theta_dot-alpha_dot
    jac = np.insert(jac, 2, jac[7, :]-jac[2, :], axis=0)                                                                     #adds a lign for gamma
    
    #New jac lign : [V, beta, gamma, alpha, p, q, r, phi, theta]
    LignLat = (1, 4, 6, 7)  #
    ColLat = (1, 4, 6, 7)  # beta, p, r, phi, delta_a
    LignLongi = (0, 2, 3, 5)  # (V, gamma, alpha, q)
    ColLongi = (0, 2, 3, 5)  # (V, gamma, alpha, q)
    TransiLat = jac[LignLat, :]
    LatJac = TransiLat[:, ColLat]
    TransiLongi = jac[LignLongi, :]
    LongiJac = TransiLongi[:, ColLongi]


    Lateigvals = scipy.linalg.eigvals(LatJac)
    Longieigvals = scipy.linalg.eigvals(LongiJac)


    print("Longitudinal Eigen value :")
    print(Longieigvals)
    print("Lateral Eigen value :")
    print(Lateigvals)






# --- display informations --- 
print("\nConfiguration : "+g.hangar['version'])
print("Vertical Tail size ratio : {0:0.2f}".format(g.VTsize))
print("Lateral stability, Cy_beta = {0:0.2f}, Cn_beta = {1:0.3f}, Cn_r = {2:0.2f}".format(CoefMatrix[1, 1], CoefMatrix[5, 1], CoefMatrix[5, 4]))
print("Default conditions : Vel_base = {0:0.1f}m/s, Beta = {1:0.1f}°, gamma={2:0.1f}0, omega={3:0.2f}°/s, Altitude={4:0.0f}m".format(V_base, beta_base/math.pi*180, gamma/math.pi*180, omega/np.pi*180,H_base) )
print("Number of engine : {0} \nNumber of inoperative engines : {1}".format(g.N_eng, g.inop))














Aero_Derivatives, Aero_Derivatives_adim = derivatives_calc.aero_coefficients(x, fixtest, Coef_base, atmospher, g, PW)



Eigenvalues_info , Eigenvalues_info_adim = Eigenvalues_manually.Eig_info(Longieigvals,Lateigvals,g,fixtest)

Eigenvalues = Eigenvalues_manually.Eigenvalues(Aero_Derivatives_adim, x, fixtest, atmospher, g, PW, CoefMatrix)




Xsample_longitudinal, CD_sample, CL_sample, Cm_sample = CL_figure.Sample_generation(x, fixtest, Coef_base, atmospher, g, PW)




start = time.time()

Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample = Apricott.Sample_generation(x, fixtest, Coef_base, atmospher, g, PW)

end = time.time()
print(end - start)


Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample = Python_to_Matlab.Python_to_Matlab(Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample)


