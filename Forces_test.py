
"""
Useful model for giving out the aerodynamic forces whichever the variables from state vector and fix vector are, so not
need for trim. If Aeroforces is changed it should be changed here too.

author: david.planas-andres

"""

import numpy as np
import math
from numpy.linalg import inv
from StabilityMapUtils import AeroForces

import ReadFileUtils as Read  # utils to read Xfoil file
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import matplotlib.pyplot as plt


def Constraints_DEP(CoefMatrix, atmo, g, PropWing):



    #x = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i(hay 12), V , beta , gamma, omega] = x + fix


    rho = atmo[1]
    PW = PropWing


    # --- Now prepare variables for equations ---
    V = 72
    beta = 0
    gamma = 0
    omega = 0


    alpha = 0.09595181246845542
    p = 0
    q = 0
    r = 0

    phi = 0
    theta = 0.09595181261150709
    aileron = 0
    elevator = -0.09167814432205841
    rudder = 0
    delta_x = 0.3496433683241823

    g.FlapDefl = 0*np.pi/180  # 15*np.pi/180 , 30*np.pi/180



    if g.FlapDefl == 0:
        g.Cd0_fl = 0
        g.CL0_fl = 0
        g.Cm0_fl = 0
        g.Cda = g.Cda_fl_0
        g.Cdb = g.Cdb_fl_0
        g.Cdc = g.Cdc_fl_0

        g.eps0 = g.eps0_flaps0
        g.deps_dalpha = g.deps_dalpha_flaps0
    elif g.FlapDefl == 15 * np.pi / 180:
        g.Cd0_fl = g.Cd0_fl_15
        g.CL0_fl = g.CL0_fl_15
        g.Cm0_fl = g.Cm0_fl_15
        g.Cda = g.Cda_fl_15
        g.Cdb = g.Cdb_fl_15
        g.Cdc = g.Cdc_fl_15

        g.eps0 = g.eps0_flaps15
        g.deps_dalpha = g.deps_dalpha_flaps15
    elif g.FlapDefl == 30 * np.pi / 180:
        g.Cd0_fl = g.Cd0_fl_30
        g.CL0_fl = g.CL0_fl_30
        g.Cm0_fl = g.Cm0_fl_30
        g.Cda = g.Cda_fl_30
        g.Cdb = g.Cdb_fl_30
        g.Cdc = g.Cdc_fl_30

        g.eps0 = g.eps0_flaps30
        g.deps_dalpha = g.deps_dalpha_flaps30


    x = np.array([alpha, p, q, r, phi, theta, aileron, elevator, rudder])

    for i in range(int(g.N_eng)):
       x = np.append(x, delta_x)

    I = np.array([[g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz]])



    # --- Compute aerodynamic forces ---
    # here subvector  must be : (alpha, beta, p, q, r, da, de,dr,  dx)
    sub_vect = np.array([alpha, beta, p, q, r])
    if g.nofin == False:
        sub_vect = np.append(sub_vect, [aileron, elevator, rudder])  # rudder is allowed
    else:
        sub_vect = np.append(sub_vect, [aileron, elevator])  # no fin allowed, default case


    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng
    Fx_vec = g.Thrust(x[-g.N_eng:],V_vect)
    Fx = np.sum(Fx_vec)



    Moment= np.zeros((g.N_eng,3))
    for i in range(g.N_eng):
        a = np.array([g.x_cg - (g.lemac - g.xp), g.PosiEng[i], g.z_m])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i + g.alpha_0+g.ip),  0,  -Fx_vec[i]*np.sin(g.alpha_i + g.alpha_0+g.ip)])
        Moment[i, :] = np.cross(a, b)
    Thrust_moment_body_axis = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))

    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta) , np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(beta)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

    Trust_moment_aero_axis = Body2Aero_matrix @  Thrust_moment_body_axis

    Mt = Trust_moment_aero_axis




    Tc = Fx_vec / (2 * rho * g.Sp * V_vect ** 2)



    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi)

    fixtest = np.array([V, beta, gamma, omega])




    if g.IsPropWing:
         h=1


    g.IsPropWing = True
    g.IsPropWingDrag = True


    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    # F contiene solo las fuerzas y momentos aerodinámicos en ejes viento.
    printx(x, fixtest, atmo, g, PW)


    A=np.zeros(10)

    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx*np.cos(alpha+g.alpha_i+g.alpha_0+g.ip)*np.cos(beta)/g.m
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m*V)
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))-Fx*np.sin(alpha+g.alpha_i+g.alpha_0+g.ip)/(g.m*V*np.cos(beta))
    A[3:6]=np.dot(inv(I), np.array([Mt[0],Mt[1],Mt[2]])+F[3:6]-np.cross(np.array([p,q,r]),np.dot(I,np.array([p,q,r]))))
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7]=q*math.cos(phi) -r*math.sin(phi)
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)




    g.IsPropWing = False
    g.IsPropWingDrag = False

    F_no_int = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)


    A=np.zeros(10)

    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx*np.cos(alpha+g.alpha_i+g.alpha_0+g.ip)*np.cos(beta)/g.m
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m*V)
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))-Fx*np.sin(alpha+g.alpha_i+g.alpha_0+g.ip)/(g.m*V*np.cos(beta))
    A[3:6]=np.dot(inv(I), np.array([Mt[0],Mt[1],Mt[2]])+F[3:6]-np.cross(np.array([p,q,r]),np.dot(I,np.array([p,q,r]))))
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7]=q*math.cos(phi) -r*math.sin(phi)
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)



    if h == 1:
        g.IsPropWing = True
        g.IsPropWingDrag = True





    return F,F_no_int

















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
        if V <= g.VelFlap:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], g.FlapDefl, g, False, beta, x[1], V, x[3])
        else:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], 0, g, False, beta, x[1], V, x[3])

    if g.nofin==False:
        print("dr = {0:0.2f}\xb0".format(x[8]/math.pi*180))






















