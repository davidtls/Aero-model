
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



    # --- Now prepare variables for equations ---
    V = 77.1677
    beta = 20*np.pi/180
    gamma = 0
    omega = 0


    alpha = 2*np.pi/180        #0.09727079748668332
    p = 0  # 0.25 max
    q = 0
    r = 0  # 0.1 max

    phi = 0
    theta = alpha
    aileron = 0
    elevator = 0
    rudder = 0
    delta_x = 1

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


    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp
    Fx_vec = g.Thrust(x[-g.N_eng:],V_vect, atmo)
    Fx = np.sum(Fx_vec)



    #Matrix to transform a vector from body reference to aero reference
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(beta)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

    #Thrust force in body reference
    F_thrust_body = [Fx*np.cos(g.alpha_i - g.alpha_0+g.ip), 0, -Fx*np.sin(g.alpha_i - g.alpha_0+g.ip)]




    # Thrust force is transformed from body to aero reference
    F_thrust_aero = Body2Aero_matrix @ F_thrust_body


    # Moment of thrust is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip), 0, -Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip)])
        Moment[i, :] = np.cross(a, b)
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))

    Mt = Thrust_moment_body



    Tc = Fx_vec / (2 * rho * g.Sp * V ** 2)





    fixtest = np.array([V, beta, gamma, omega])

    sinbank = np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank = np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi)


    if g.IsPropWing:
         h = 1


    g.IsPropWing = False
    g.IsPropWingDrag = False


    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    # F contiene solo las fuerzas en ejes viento y momentos aerodinámicos en ejes cuerpo.

    printx(x, fixtest, atmo, g, PropWing)


    CL = -F[2]/(0.5*rho*V**2 * g.S)
    CD = -F[0]/(0.5*rho*V**2 * g.S)
    CY = -F[1]/(0.5*rho*V**2 * g.S)
    Clroll = -F[3]/(0.5*rho*V**2 * g.S * g.b)
    Cm = -F[4]/(0.5*rho*V**2 * g.S * g.c)
    Cn = -F[5]/(0.5*rho*V**2 * g.S * g.b)


    A=np.zeros(10)

    A[0] = -9.81*np.sin(gamma)+(F[0]+F_thrust_aero[0])/g.m
    A[1] = (p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + (F[1]+ F_thrust_aero[1])/(g.m*V)
    A[2] = -(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta) + 9.81*cosbank/(V*np.cos(beta)) + (F[2] + F_thrust_aero[2])/(g.m*V*np.cos(beta))
    A[3:6] = np.dot(inv(I), np.array([Mt[0], Mt[1], Mt[2]])+F[3:6]-np.cross(np.array([p, q, r]), np.dot(I, np.array([p, q, r]))))
    A[6] = p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7] = q*math.cos(phi) - r * math.sin(phi)
    A[8] = -np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9] = -omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)




    g.IsPropWing = False
    g.IsPropWingDrag = False

    F_no_int = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)


    A=np.zeros(10)

    A[0] = -9.81*np.sin(gamma)+(F[0]+F_thrust_aero[0])/g.m
    A[1] = (p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + (F[1]+ F_thrust_aero[1])/(g.m*V)
    A[2] = -(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta) + 9.81*cosbank/(V*np.cos(beta)) + (F[2] + F_thrust_aero[2])/(g.m*V*np.cos(beta))
    A[3:6] = np.dot(inv(I), np.array([Mt[0], Mt[1], Mt[2]])+F[3:6]-np.cross(np.array([p, q, r]), np.dot(I, np.array([p, q, r]))))
    A[6] = p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7] = q*math.cos(phi) - r * math.sin(phi)
    A[8] = -np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9] = -omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)



    if h == 1:
        g.IsPropWing = True
        g.IsPropWingDrag = True





    return F,F_no_int






def Constraints_DEP_body(CoefMatrix, atmo, g, PropWing):
    """function defining constraints for power minimization
    inputs:
    -x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
    x is the state to determine
    length of x except the propulsion levels is 8
    -fix = [V, beta, gamma, omega]
    fix is the vector of parameters whom are fixed by the user

    """


    #x = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i(hay 12), V , beta , gamma, omega] = x + fix


    rho = atmo[1]
    PW = PropWing


    # --- Now prepare variables for equations ---
    V = 72
    beta = 0
    gamma = 0
    omega = 0


    alpha = 0.09727079748668332
    p = -3.5734202462290796e-22  # 0.25 max
    q = 8.205631676526035e-22
    r = -3.6545873438358454e-22   # 0.1 max

    phi = -9.355764297366552e-06
    theta = 0.09727079748245356
    aileron = -8.412261477902928e-05
    elevator = -0.09318168652099662
    rudder = -2.3758828310177916e-06
    delta_x = 0.3153640459583273

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
    #here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect = np.array([alpha, beta, p, q, r])
    if g.nofin == False:
        sub_vect = np.append(sub_vect, [x[6], x[7], x[8]])  # rudder is allowed
    else:
        sub_vect = np.append(sub_vect, [x[6], x[7]])  # no fin allowed, default case



    #Thrust forces and moments

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp


    Fx_vec = g.Thrust(x[-g.N_eng:], V_vect, atmo)
    Fx = np.sum(Fx_vec)


    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)                                                                                       #For adimension V, has already been used for calculating FXi

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    #F gives out aerodinamical forces in aero axis: Drag, lateral force and lift and moments
    # Does not give out X,Y,Z





    #Matrix to transform a vector from body reference to aero reference
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(beta)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

    #Thrust force in body reference
    F_thrust_body = [Fx*np.cos(g.alpha_i - g.alpha_0+g.ip) , 0 , -Fx*np.sin(g.alpha_i - g.alpha_0+g.ip)]



    # Moment of thrust is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip), 0,-Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip)])
        Moment[i, :] = np.cross(a, b)
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))

    Mt = Thrust_moment_body



    # Transformation of aerodynamic forces from aero to body reference
    F_aero_body=np.zeros(int(len(F)))
    F_aero_body[0:3] = np.transpose(Body2Aero_matrix) @ F[0:3]



    # Transformation of aerodynamic speed from aero to body reference
    [u,v,w] = np.transpose(Body2Aero_matrix) @ np.concatenate(([V], [0], [0]))


    sinbank = np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank = np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi)


    A = np.zeros(10+g.inop)

    A[0] = (1/g.m) * (F_thrust_body[0] + F_aero_body[0]) - 9.81*np.sin(theta) +r*v - q*w
    A[1] =(1/g.m)*(F_aero_body[1]) + 9.81*np.cos(theta)*np.sin(phi) - r*u + p*w
    A[2] =(1/g.m)*(F_thrust_body[2] +F_aero_body[2])+ 9.81*np.cos(theta)*np.cos(phi) +q*u - p*v

    A[3:6] = np.dot(inv(I), np.array([Mt[0], Mt[1], Mt[2]])+F[3:6]-np.cross(np.array([p, q, r]), np.dot(I, np.array([p, q, r]))))

    A[6] = p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7] = q*math.cos(phi) - r * math.sin(phi)

    A[8] = -np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9] = -omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)



    for i in range(g.inop):
        A[-1-i] = x[-1-i]                                                                                                 #The inoperative engine are the last ones (right wing). Its value is minimized (to zero)

    if g.hangar['version'] == 'original':                                                                                 #For obligating all the engines to have the same thrust
        #no DEP with original twin or N engines; all engines have the same thrust
        D = np.copy(A)
        for i in range(g.N_eng-g.inop-1):
            AAd = x[-g.N_eng]-x[-g.N_eng+i+1]
            D = np.append(D, [AAd])
        return D
    else:
        return A

















def Long_equilibrium2(CoefMatrix, atmo, g, PropWing):
    """function defining constraints for speed minimization in longitudinal
    inputs:
        -x =[V, alpha, theta, delta_e, delta_i]
        x is the state to determine
        length of x except the propulsion levels is 8
        -fix = [gamma, beta, p, q, r, phi, da, dr]
        fix is the vector of parameters whom are fixed by the user

        Again
        gamma = beta = p = q = r = phi = da = dr = 0 as we are in LONGITUDINAL equilibrium

        4 equations (2 forces, 1 moment, theta = alpha + gamma)
        (variables = V, alpha, theta, de, di) problem oversized
         that means there is place for optimization, with objective function V)

    """


    rho = atmo[1]



    # --- Now prepare variables for equations ---
    V =  70
    alpha = 0.10723828379856166
    de =  -0.10514575695526691
    dx =   0.30905799124088



    beta = 0
    gamma = 0

    p = 0
    q = 0
    r = 0

    phi = 0
    theta = alpha
    da = 0

    dr = 0


    g.FlapDefl = 0 * np.pi/180  # 15*np.pi/180 , 30*np.pi/180



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




    I = np.array([[g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz]])


    # --- Compute aerodynamic forces ---
    #here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])  # rudder is allowed


    #Thrust forces and moments

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp



    Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
    Fx = np.sum(Fx_vec)



    #Matrix to transform a vector from body reference to aero reference
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(beta)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

    #Thrust force in body reference
    F_thrust_body = [Fx*np.cos(g.alpha_i - g.alpha_0+g.ip), 0, -Fx*np.sin(g.alpha_i - g.alpha_0+g.ip)]




    # Thrust force is transformed from body to aero reference
    F_thrust_aero = Body2Aero_matrix @ F_thrust_body


    # Moment of thrust is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip), 0,-Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip)])
        Moment[i, :] = np.cross(a, b)
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))

    Mt = Thrust_moment_body


    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)                                                                                       #For adimension V, has already been used for calculating FXi

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)


    #F gives out aerodinamical forces in aero axis: Drag, lateral force and lift and moments
    # Does not give out X,Y,Z


    #     Now sum up the constraints:


    A = np.zeros(4)


    A[0] = +(F[0] + F_thrust_aero[0])
    A[1] = 9.81*g.m + (F[2] + F_thrust_aero[2])
    A[2] = (Mt[1] + F[4])
    A[3] = alpha + gamma - theta


    fixtest = np.array([V, beta, gamma, 0])

    x = np.array([alpha, p, q, r, phi, theta, da, de, dr])

    for i in range(int(g.N_eng)):
         x = np.append(x, dx)

    #printx(x, fixtest, atmo, g, PropWing)

    CL = -F[2]/(0.5*rho*V**2 * g.S)
    CD = -F[0]/(0.5*rho*V**2 * g.S)





    return A
























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

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * fix[1] + g.wingsweep) - x[3] * g.yp

    if g.IsPropWing:
        if V <= g.VelFlap or g.FlapDefl != 0:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], g.FlapDefl, g, False, beta, x[1], V, x[3])
        else:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], 0, g, False, beta, x[1], V, x[3])

    if g.nofin==False:
        print("dr = {0:0.2f}\xb0".format(x[8]/math.pi*180))






















