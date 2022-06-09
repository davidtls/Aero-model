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
import scipy.io #input/output with matlab
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import fsolve
import sys
sys.path.insert(0, '/home/e.nguyen-van/Documents/codesign-small-tail/Python/PattersonAugmented')
sys.path.insert(0, '/home/e.nguyen-van/Documents/codesign-small-tail/PythonStabilityMapUtils')
sys.path.insert(0, '/home/e.nguyen-van/Documents/codesign-small-tail/Python/AircraftClass')
sys.path.insert(0, '/home/e.nguyen-van/Documents/codesign-small-tail/Python')



import PattersonAugmented as PA
import control
import pylab
from StabilityMapUtils import ReadFileUtils
from StabilityMapUtils import equation as e
from StabilityMapUtils import AeroForces
import time
import pickle
import datetime
import Eigenvalues_manually
import derivatives_calc
import Apricott
import Python_to_Matlab


import Forces_test
#from StabilityMapUtils import equation_body as e    # cant be uncommented line 33 is uncommented


"""
Things to check:
aircraft model
Propeller Wing Interaction
inop engines
Speed
Altitude
Complete Optimization/ Long_equilibrium/ Long_Equilibrium2
Forces comparison deactivated
g.hangar DEPoriginal (different thrust) or original (same thrust)
flaps
"""






aircraft = {'model': 'ATR'}   #OPTIONS: ATR, DECOL






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
    V_base = 72          #1.3 * Vsr
    beta_base = 0 / 180 * math.pi
    gamma = 0                                                                                                           # previous condition  np.arctan(3 / 100)  # (3/100)/180*np.pi##math.atan(0/87.4)#/180*math.pi # 3% slope gradient # Best climb rate: 6.88m/s vertical @ 87.5m/s = 4.5°gamma, see http://www.atraircraft.com/products_app/media/pdf/Fiche_72-600_Juin-2014.pdf
    R = 000  # in meters the turn radius
    phimax = 5  # in degree the max bank angle authorized
    alphamax = 25  # in degree, stall bound for trimming
    alphastall = 11.7  # that's for patterson
    deltaRmax = 30  # in degree
    ThrottleMax = 1  # max thrust level
    ThrottleMin = 1e-9  # min throttle, don't accept 0 thrust
    ThrottleMinExt = 1e-9  # -0.34


    # --- dictionnary for type of aircraft studied. aircraft: ATR72, version : 'original', 'DEPoriginal', 'DEPnofin'
    g.hangar = {'aircraft': 'ATR72', 'version': 'original'}

    g.VelFlap = 0 # modificada para que no salte   71  # in m/s the maximum velocity at which flap are deployed
    g.CL0_fl = g.CL0_fl  # / 2   # for take off in no-interaction divide by two








    if g.IsPropWing:
        # ensures alpha fuselage as reference for stall
        g.alpha_max = alphastall / 180 * np.pi + g.alpha_i - g.ip
        g.alpha_max_fl = alphastall / 180 * np.pi + g.alpha_i - g.alpha_0
    else:
        g.alpha_max = alphastall / 180 * np.pi
        g.alpha_max_fl = alphastall / 180 * np.pi

    # --- additional parameters (default edited during execution) ---
    g.set_nofin(False)  # =    True means : no rudder used    False:rudder used




    """ Algorithm set up """
    # ---- Optim parameter ------
    MaxIter = 100
    tolerance = 1e-5
    method = 'SLSQP'  # 'trust-interior' or 'SLSQP'

    g.Pkeyword = 'DefaultPatterson'










mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif'] = ['Latin Modern Roman'] #for legend / FreeSerif
mpl.rcParams["mathtext.fontset"] = "stix" # for math in legend
mpl.rcParams["font.size"] = 16 






print("Engine position :")
strout = ""
for i in range(len(g.PosiEng)):
    strout = strout+str(i+1)+" : "+"{:.3}".format(g.PosiEng[i])+"m, "
print(strout)





g.Complete_Optimization = True
g.Long_equilibrium = False
g.Long_equilibrium2 = False




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

       Velocities = (70, 90 , 110  , 130 , 150)                                                                          #Two first speeds for H = 0 m, last 3 for H=5000m
       rho_vec = (1.225 , 1.225, 0.736116,  0.736116,  0.736116)                                               # rho= 0.736116 kg/m^3     H=5000 m   a=320.529 m/s
       #a_sound=(340,340,320.529 ,320.529 ,320.529)
       Mach = [0.2058, 0.2647, 0.3431, 0.4055, 0.4679]                                                         # Mach= Vel/a




elif aircraft['model']=='DECOL':

      Velocities = (10, 15, 20, 25, 30, 35)
      rho_vec = (1.225, 1.225, 1.225, 1.225, 1.225, 1.225)
      Mach = np.ones((len(Velocities), 1))*0.0001












""" Algorithm start """

#--- List all .stab file from vsp aero and read the coeff ---- 

if g.hangar['aircraft'] == 'ATR72':
    if g.hangar['version'] == 'original' or g.hangar['version']=='DEPoriginal':
        if g.hangar['version'] == 'original' and g.nofin == True:
            print("WARNING : Using "+g.hangar['version']+" without rudder. Not recommended...")
        if g.hangar['version'] == 'DEPoriginal' and g.inop != 0 and g.nofin == True:
            print("WARNING : Using "+g.hangar['version']+" with inoperative engine and without rudder. Not recommended...")
#        filename=['ATR72_SI_MTOW_Control_flap.stab', 'ATR72_SI_MTOW_Control_flapVmax.stab',
#                  'ATR72_SI_MTOW_Control_Vmin.stab','ATR72_SI_MTOW_Control_Manoeuver.stab',
#                  'ATR72_SI_MTOW_Control_Cruise.stab']
        path = 'ATR72_SI_MTOW_FinLess_STAB/'
        filenameNoFin = [path+'ATR72_FinLess_mach1.stab', path+'ATR72_FinLess_mach2.stab', path+'ATR72_FinLess_mach3.stab',
                   path+'ATR72_FinLess_mach4.stab', path+'ATR72_FinLess_mach5.stab']
        Matrix = ReadFileUtils.ReadStabCoef(filenameNoFin)
        

elif g.hangar['aircraft'] == 'DECOL':

          # --- List all .stab file from vsp aero and read the coeff ----
          path = 'DECOL_STAB/'  # 'home/e.nguyen-van/Documents/DECOLStability&Analysis/DECOLDATA/DECOLGeom_DegenGeom_6_3_18h14
          filenameNoFin = [path + '_FinLess_Vinf10000.stab',
                           path + '_FinLess_Vinf15000.stab',
                           path + '_FinLess_Vinf20000.stab',
                           path + '_FinLess_Vinf25000.stab',
                           path + '_FinLess_Vinf30000.stab',
                           path + '_FinLess_Vinf35000.stab']
          MatrixNoFin = ReadFileUtils.ReadStabCoef(filenameNoFin)

          # copy the matrix to avoid error and keep track
          Matrix = np.copy(MatrixNoFin)




print("Adjusting Kf, Kh and VT size")
print("New Kf and Kh")
print(g.AdjustVT())

CoefMatrix=g.NicolosiCoef(Matrix[:, 1:], Mach)
        

# --- Interpol coefficients for test case ---
# Find sound velocity and air density
atmospher = g.GetAtmo(H_base)
a_sound = atmospher[0]
rho_base = atmospher[1]
M_base = V_base/a_sound







if g.hangar['aircraft'] == 'ATR72':
               Coef_base = AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)
               g.Matrix_no_tail_terms = AeroForces.CoefInterpol(M_base, Matrix[:, 1:], Mach)

elif g.hangar['aircraft'] == 'DECOL':
               Coef_base = AeroForces.CoefInterpol(V_base, CoefMatrix, Velocities)
               g.Matrix_no_tail_terms = AeroForces.CoefInterpol(M_base, Matrix[:, 1:], Velocities)







# Define here the PropWing interaction

if aircraft['model'] =='ATR':

        PropPath = "./ATR72_SI_MTOW_Control_FinLess_FEM/"
        PropFilenames = {'fem':[PropPath+"ATR72_FinLess_mach1",
                                PropPath+"ATR72_FinLess_mach2",
                                PropPath+"ATR72_FinLess_mach3",
                                PropPath+"ATR72_FinLess_mach4",
                                PropPath+"ATR72_FinLess_mach5"],
                        'AirfoilPolar':PropPath+"naca3318Pol.txt",
                        'FlapPolar':PropPath+"naca3318fl+10.txt",
                        'AileronPolar':PropPath+"naca3318fl+10.txt"} # format for prop file : [[Cldist=f(M)],polar clean airfoil, polar flap, polar aile]
        g.PolarAilDeflDeg = 10
        g.PolarFlDeflDeg = 10
        g.AilDiff = 0.5

        if g.IsPropWing:
                     g.Pkeyword = "DefaultPatterson"
        PW = PA.PropWing(g, PropFilenames)
        PW.DeltaCL_a_0 = 1




elif aircraft['model']=='DECOL':

            g.PolarFlDeflDeg = 5
            g.PolarAilDeflDeg = 5
            PropPath = "DECOL_FEM/"
            PropFilenames = {'fem':[PropPath+"_FinLess_Vinf10000.0"],
                            'AirfoilPolar':PropPath+"S3010_XTr10_Re350.txt",
                            'FlapPolar':PropPath+"S3010_XTr10_Re350_fl5.txt",
                            'AileronPolar':PropPath+"S3010_XTr10_Re350_fl5.txt"} # format for prop file : [[Cldist=f(M)],polar clean airfoil, polar flap, polar aile]
            PW = PA.PropWing(g, PropFilenames)
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
#Forces_comparison = Forces_test.Long_equilibrium2(Coef_base, atmospher, g, PW)





if g.Complete_Optimization == True:


        #Initialise test and guess vectors
        if g.nofin == False:
            # x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
            x0 = np.array([5*math.pi/180, 0, 0, 0, 0.00, 0.0, 0.0, 0.0, 0])
            bnds = ((-5*math.pi/180, alphamax*math.pi/180), (-0.2, 0.2), (-0.2, 0.2), (-0.2, 0.2), (-phimax/180*math.pi, phimax/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-20/180*math.pi, 20/180*math.pi), (-deltaRmax/180*math.pi, deltaRmax/180*math.pi))
        else:
            # x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_i]
            x0 = np.array([7.5*math.pi/180, 0, 0, 0, 0.00, 0.05, 0.00, 0.00])
            bnds = ((-5*math.pi/180, alphamax*math.pi/180), (-0.2, 0.2), (-0.2, 0.2), (-0.2, 0.2), (-phimax/180*math.pi, phimax/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-30/180*math.pi, 30/180*math.pi), (-20/180*math.pi, 20/180*math.pi))

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









if g.Long_equilibrium == True:

        g.IsPropWing = True
        g.IsPropWingDrag = True
        #Initialise test and guess vectors
        # x =[V, theta, delta_e]
        x0 = np.array([60, 0, 0])




        alpha = 8 * np.pi/180
        dx = 0.5
        dx_vec = np.full(g.N_eng, dx)

        beta = 0
        p = 0
        q = 0
        r = 0
        phi = 0
        da = 0
        dr = 0
        omega = 0


        fixtest = np.array((alpha,p,q,r,beta,phi,da,dr))

        fixtest=np.concatenate((fixtest, dx_vec))

        # put everything in tuples for passing to functions
        diccons = (np.copy(fixtest), np.copy(Coef_base), atmospher, g, PW)

        # --- solving algorithm ---

        k = fsolve(e.Long_equilibrium, x0, args=diccons, xtol=1e-10, maxfev=0, factor=0.1)


        # check if constraints are validated
        constraints_calc = e.Long_equilibrium(k, *diccons)
        print("\nConstraints")


        #Reorders vector for leaving them as in the minimization in power case

        gamma = k[1] - alpha
        alpha = alpha
        V = k[0]


        #Reorganize everything
        x1 = np.array((alpha, p, q, r, phi, k[1], da, k[2], dr))
        x = np.concatenate((x1, fixtest[-g.N_eng:]))
        fixtest = np.array([V, beta, gamma, omega])
        diccons = (np.copy(fixtest), np.copy(Coef_base), atmospher, g, PW)

        # for coherance shall be x = [alpha, p, q ,r, phi, theta, da, de, dr, dx]
        #                  fixtest = [V, beta, gamma, omega]



if g.Long_equilibrium2 == True:




    #Initialise test and guess vectors
    # x =[V, alpha, theta, delta_e, delta_x]
    x0 = np.array([20, 5*math.pi/180, 0, 0, 0.4])

    if g.IsPropWing:
         bnds = ((20, 100), (-10*math.pi/180, 25*math.pi/180), (-10/180*math.pi, 25/180*math.pi), (-23/180*math.pi, 13/180*math.pi), ((ThrottleMin, ThrottleMax)))
    else:
         bnds = ((20, 100), (-10*math.pi/180, 0.2561481435699163), (-10/180*math.pi, 25/180*math.pi), (-23/180*math.pi, 13/180*math.pi), ((ThrottleMin, ThrottleMax)))

    #For maximizing speed (instead of minimizing it for stall speed), function turn to 1/f, anyway is not well defined so
    #    x0 = np.array([150, 5*math.pi/180, 0, 0, 0.99]) and
    #   bnds = ((100, 200), (-10*math.pi/180, 25*math.pi/180), (-10/180*math.pi, 25/180*math.pi), (-23/180*math.pi, 13/180*math.pi), ((0.96, ThrottleMax)))


    # --- imposed conditions ---
    gamma = 0
    beta = 0
    p = 0
    q = 0
    r = 0
    phi = 0
    da = 0
    dr = 0
    omega = 0


    fixtest = np.array([gamma, beta, p, q, r, phi, da, dr, omega])


    # put everything in tuples for passing to functions
    diccons = (np.copy(fixtest), np.copy(Coef_base), atmospher, g, PW)  # fix, CoefMatrix,Velocities, rho, g




    # --- minimization algorithm ---
    ## SLSQP
    if method == 'SLSQP':
        k = minimize(e.V_min, np.copy(x0), args=diccons, bounds=bnds,
                     constraints={'type': 'eq', 'fun': e.Long_equilibrium2, 'args': diccons},
                     options={'maxiter': MaxIter, 'disp': True}, tol=tolerance)
    elif method == 'trust-interior':
        ## Trust interior
        k = minimize(e.V_min, np.copy(x0), args=diccons, method='trust-constr', bounds=bnds,
                     constraints={'type': 'eq', 'fun': e.Long_equilibrium2, 'args': diccons},
                     options={'maxiter': MaxIter, 'disp': True}, tol=tolerance)
    else:
        sys.exit('Non valid optimization method')

    #leaving things like before

    constraints_calc = e.Long_equilibrium2(k.x, *diccons)
    print("\nConstraints")


    x1 = np.array((k.x[1], p, q, r, phi, k.x[2], da, k.x[3], dr))
    x = np.concatenate((x1, k.x[-g.N_eng:]))

    fixtest = np.array([k.x[0], beta, gamma, omega])
    diccons = (np.copy(fixtest), np.copy(Coef_base), atmospher, g, PW)

    # for coherance shall be x = [alpha, p, q ,r, phi, theta, da, de, dr, dx]
    #                  fixtest = [V, beta, gamma, omega]








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
    print("da={0:0.2f}\xb0, de= {1:0.2f}\xb0".format(da,de))

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * fix[1] + g.wingsweep) - x[3] * g.PosiEng

    if g.IsPropWing:
        if V <= g.VelFlap or g.FlapDefl != 0:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], g.FlapDefl, g, False, fix[1], x[1], V, x[3])
        else:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], 0, g, False, fix[1], x[1], V, x[3])

    if g.nofin==False:
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




start = time.time()

Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample = Apricott.Sample_generation(x, fixtest, Coef_base, atmospher, g, PW)

end = time.time()
print(end - start)


Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample = Python_to_Matlab.Python_to_Matlab(Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample)


