
"""
Useful model for computing the aerodynamic forces whichever the variables from state vector and fix vector are.
Therefore they can be calculated outside a trim state.
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




def Distributed_engines(N_eng, DiameterHLP ):
    """
    Function valid for the X-57 only. This function receives a number of engines and their diameter and distributes them equally
    among the available span that starts outside the fuselage and finish at a distance of the wingtip equal to
    1 radius of the wing-tip cruise propeller. It returns the value of yp, the vector containing the y position of
    the engines (from left wing to right wing) in the body
    reference system.

    Modify so that for more engines the diameter is reduced.
    """

    b = 9.642  # wingspan
    Rtip = 0.5 * 60 * 0.0254  # radius of propeller at tip
    FusR = 0.60198  # max. radius fuselage

    change_diameter = True

    if N_eng <= 12:

         if change_diameter:
             """
             For less than 12 engines, we increase the diameter in order to have still or the wing area washed.
             """
             avail_y = 0.5*b - Rtip - FusR
             RHLP = avail_y/(2*N_eng/2)

             max_y = (0.5*b - Rtip - RHLP)
             min_y = (FusR + RHLP)

             yp = np.linspace(min_y, max_y, int(N_eng/2))
             yp = np.hstack((np.flip(-yp), yp))


         else:
             # Diameter is not modified

             RHLP = 0.5*DiameterHLP       # radius of HLP

             max_y = (0.5*b - Rtip - RHLP)
             min_y = (FusR + RHLP)

             yp = np.linspace(min_y, max_y, int(N_eng/2))
             yp = np.hstack((np.flip(-yp), yp))

    else:  # Diameter needs to be modified so that there is space. In this case is modified for keeping all area wet

        avail_y = 0.5*b - Rtip - FusR
        RHLP = avail_y/(2*N_eng/2)

        max_y = (0.5*b - Rtip - RHLP)
        min_y = (FusR + RHLP)

        yp = np.linspace(min_y, max_y, int(N_eng/2))
        yp = np.hstack((np.flip(-yp), yp))

    return yp, RHLP





def Constraints_DEP(CoefMatrix, atmo, g, PropWing):



    #x = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i(hay 12), V , beta , gamma, omega] = x + fix


    rho = atmo[1]



    # --- Now prepare variables for equations ---
    V = 70  # 51.741863715877166
    beta = 0
    gamma = 0
    omega = 0


    alpha = 0.2362537139286036
    p = 0  # 0.25 max
    q = 0
    r = 0  # 0.1 max

    phi = 0
    theta = alpha
    aileron = 0
    elevator = -0.20600870692395656
    rudder = 0
    delta_x = 0.4    # dx = 0.5


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

    #engines = np.array([0.5,0.5,0.5,0.8,0.5,0.5,0.5,0.5,0.2,0.5,0.5,0.5])
    #x = np.append(x, engines)

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
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])


    # Moment and Force of thrust  is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    F_thrust_body = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip[i]), 0,-Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip[i])])
        Moment[i, :] = np.cross(a, b)
        F_thrust_body[i,:] = b
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))
    F_thrust_body = np.array((np.sum(F_thrust_body[:, 0]), np.sum(F_thrust_body[:, 1]), np.sum(F_thrust_body[:, 2])))

    Mt = Thrust_moment_body

    # Thrust force is transformed from body to aero reference
    F_thrust_aero = Body2Aero_matrix @ F_thrust_body




    Tc = Fx_vec / (2 * rho * g.Sp * V ** 2)





    fixtest = np.array([V, beta, gamma, omega])

    sinbank = np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank = np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi)


    if g.IsPropWing:
         h = 1


    g.IsPropWing = True
    g.IsPropWingDrag = True


    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    # F contiene solo las fuerzas en ejes viento y momentos aerodinámicos en ejes cuerpo.

    printx(x, fixtest, atmo, g, PropWing)


    #CL = -F[2]/(0.5*rho*V**2 * g.S)
    #CD = -F[0]/(0.5*rho*V**2 * g.S)
    #CY = -F[1]/(0.5*rho*V**2 * g.S)
    #Cm = -F[4]/(0.5*rho*V**2 * g.S * g.c)






    Croll_thrust = Mt[0]/ (0.5 * rho * V ** 2 * g.S * g.b)
    Cyawthrust = Mt[2]/ (0.5 * rho * V ** 2 * g.S * g.b)

    Clroll = F[3]/(0.5*rho*V**2 * g.S * g.b)
    Cn = F[5]/(0.5*rho*V**2 * g.S * g.b)




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
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])


    # Moment and Force of thrust  is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    F_thrust_body = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip[i]), 0,-Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip[i])])
        Moment[i, :] = np.cross(a, b)
        F_thrust_body[i,:] = b
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))
    F_thrust_body = np.array((np.sum(F_thrust_body[:, 0]), np.sum(F_thrust_body[:, 1]), np.sum(F_thrust_body[:, 2])))

    Mt = Thrust_moment_body

    # Thrust force is transformed from body to aero reference
    F_thrust_aero = Body2Aero_matrix @ F_thrust_body


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

    """
    X-57 equilibrium conditions
    V = 40.145
    alpha = 0.0070674911658034625
    de = -0.31381389162636214
    dx = 0.9999996478122892    
    """

    rho = atmo[1]

    # --- Now prepare variables for equations ---
    V = 40.125
    alpha = 8*np.pi/180
    de = 0
    dx = 0  # for zero thrust 0.022628 only at v=40.125


    beta = (0*np.pi/180)
    gamma = 0
    p = (0*np.pi/180)
    q = 0
    r = (0*np.pi/180)
    phi = 0
    theta = alpha
    da = 0
    dr = 0
    g.FlapDefl = 30 * np.pi/180  # 15*np.pi/180 , 30*np.pi/180


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

    elif g.FlapDefl == 10 * np.pi / 180:
        g.Cd0_fl = g.Cd0_fl_30
        g.CL0_fl = g.CL0_fl_30
        g.Cm0_fl = g.Cm0_fl_30
        g.Cda = g.Cda_fl_30
        g.Cdb = g.Cdb_fl_30
        g.Cdc = g.Cdc_fl_30

        g.eps0 = g.eps0_flaps30
        g.deps_dalpha = g.deps_dalpha_flaps30

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
    # here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])  # rudder is allowed


    # Thrust forces and moments

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp



    Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
    Fx = np.sum(Fx_vec)





    # Matrix to transform a vector from body reference to aero reference
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

    # Moment and Force of thrust  is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    F_thrust_body = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip[i]), 0, -Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip[i])])
        Moment[i, :] = np.cross(a, b)
        F_thrust_body[i, :] = b
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))
    F_thrust_body = np.array((np.sum(F_thrust_body[:, 0]), np.sum(F_thrust_body[:, 1]), np.sum(F_thrust_body[:, 2])))

    Mt = Thrust_moment_body

    # Thrust force is transformed from body to aero reference
    F_thrust_aero = Body2Aero_matrix @ F_thrust_body





    # Convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)

    # F gives CD CY CL aerodynamic forces in aero ref, Moments in body ref.

    CD = -F[0]/(0.5*rho*V**2 * g.S)
    CY = -F[1]/(0.5*rho*V**2 * g.S)
    CL = -F[2]/(0.5*rho*V**2 * g.S)
    Clroll = F[3]/(0.5*rho*V**2 * g.S * g.b)
    Cm = F[4]/(0.5*rho*V**2 * g.S * g.c)
    Cn = F[5]/(0.5*rho*V**2 * g.S * g.b)

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

    printx(x, fixtest, atmo, g, PropWing)









    # Plotting
    alpha_vector = np.linspace(-4, 24, 15)*np.pi/180
    dx_vector = np.linspace(0.39, 0.39, 1)
    Tc = np.zeros((len(dx_vector), g.N_eng))

    Cm_matrix = np.zeros((len(dx_vector), len(alpha_vector)))
    CL_matrix = np.zeros((len(dx_vector), len(alpha_vector)))
    CD_matrix = np.zeros((len(dx_vector), len(alpha_vector)))
    epsilon_matrix = np.zeros((len(dx_vector), len(alpha_vector)))


    Cdi = np.zeros((len(dx_vector), len(alpha_vector)))
    tempCdo = np.zeros((len(dx_vector), len(alpha_vector)))
    Cdwash = np.zeros((len(dx_vector), len(alpha_vector)))


    for i in range(len(dx_vector)):
        Fx_vec = g.Thrust(np.full(g.N_eng, dx_vector[i]), V_vect, atmo)
        Tc[i, :] = Fx_vec/(2*rho*g.Sp*V**2)

        for j in range(len(alpha_vector)):
            sub_vect = np.array([alpha_vector[j], beta, p, q, r, da, de, dr])
            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc[i, :], atmo, g, PropWing)

            Drags = PropWing.CalcCoef(Tc[i, :], V/atmo[0], atmo, alpha_vector[j], 0, g.FlapDefl, g, beta, p, V, r)

            Cdi[i, j] = Drags[2]
            tempCdo[i, j] = Drags[3]
            Cdwash[i, j] = Drags[5]

            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)
            CD_matrix[i, j] = -F[0] / (0.5*rho*V**2 * g.S)
            Cm_matrix[i, j] = F[4] / (0.5*rho*V**2 * g.S * g.c)
            epsilon_matrix[i, j] = AeroForces.Cm_and_CL_tail(V, CoefMatrix, sub_vect, Tc[i, :], atmo, g, PropWing)[2]

    print('alpha')
    print(alpha)
    print('CL')
    print(CL_matrix)
    print('Cdi')
    print(Cdi)
    print('tempCd0')
    print(tempCdo)
    print('Cdwash')
    print(Cdwash)
    print('CdTOTAL')
    print(CD_matrix)
    print('Cm')
    print(Cm_matrix)

    fig1 = plt.figure()
    ax1 = fig1.gca()
    for i in range(len(dx_vector)):
       ax1.plot(alpha_vector*180/np.pi, Cm_matrix[i, :], label="$T_c$ = {0:0.3f}".format(Tc[i, 0]), linestyle=":", color='r', alpha = 0.8 + 0.8*i/11)

    ax1.set_xlabel('alpha (°)')
    ax1.set_ylabel('Cm')
    ax1.legend()
    ax1.grid()
    fig1.tight_layout()


    fig2 = plt.figure()
    ax2 = fig2.gca()
    for i in range(len(dx_vector)):
        ax2.plot(alpha_vector*180/np.pi, epsilon_matrix[i, :]*180/np.pi, label="$T_c$ = {0:0.3f}".format(Tc[i, 0]), linestyle=":", color='g', alpha = 0.8 + 0.8*i/11)

    ax2.set_xlabel('alpha (°)')
    ax2.set_ylabel('Downwash')
    ax2.legend()
    ax2.grid()
    fig2.tight_layout()


    fig3 = plt.figure()
    ax3 = fig3.gca()
    for i in range(len(dx_vector)):
        ax3.plot(alpha_vector*180/np.pi, CL_matrix[i, :], label="$T_c$ = {0:0.3f}".format(Tc[i, 0]), linestyle=":", color='g', alpha = 0.8 + 0.8*i/6)

    ax3.set_xlabel('alpha (°)')
    ax3.set_ylabel('CL')
    ax3.legend()
    ax3.grid()
    fig3.tight_layout()



    fig4 = plt.figure()
    ax4 = fig4.gca()
    for i in range(len(dx_vector)):
        ax4.plot(alpha_vector*180/np.pi, CD_matrix[i, :], label="$T_c$ = {0:0.3f}".format(Tc[i, 0]), linestyle=":", color='b', alpha = 0.8 + 0.8*i/6)

    ax4.set_xlabel('alpha (°)')
    ax4.set_ylabel('CD')
    ax4.legend()
    ax4.grid()
    fig4.tight_layout()

    plt.show(block=True)





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
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect, atmo)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], g.FlapDefl, g, False, beta, x[1], V, x[3])
        else:
            PW.PlotDist(g.Thrust(x[-g.N_eng:], V_vect, atmo)/(2*atmo[1]*g.Sp*V**2), V/atmo[0], atmo, x[0], x[6], 0, g, False, beta, x[1], V, x[3])

    if g.nofin==False:
        print("dr = {0:0.2f}\xb0".format(x[8]/math.pi*180))



























def PLOTSX57(CoefMatrix, atmo, g, PropWing):

    """
    % Aircraft: X-57
    % H = 0 m, V = 28.3 m/s , beta = 0, gamma = 0, omega = 0, dx = 1
    % Flaps deployed: n
    % dx = 1
    % de = 0 (But we will just plot CL,wing anyways)

    """

    rho = atmo[1]
    # --- Now prepare variables for equations ---
    V = 28.3 #40.145
    de = 0
    dx = 1
    beta = 0
    gamma = 0
    p = 0
    q = 0
    r = 0
    phi = 0
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

    ang_speed = np.array([p, q, r])


    # Thrust forces and moments
    #V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp
    original_ip_vector = g.ip
    original_offset_vector = g.x_offset
    original_diameter = g.Dp


    #PLOT 1
    #Angles of attack: 4, 8 [°] (two lines)
    #Axe x: Variation of installation angle
    #Axe y: CL,wing

    """
    alpha_vector = np.linspace(8, 8, 1)*np.pi/180
    ip_vector = np.linspace(g.ip[0]*180/np.pi-20, g.ip[0]*180/np.pi+8, 29)*np.pi/180
    CL_matrix = np.zeros((len(alpha_vector), len(ip_vector)))

    for i in range(len(alpha_vector)):

        for j in range(len(ip_vector)):
            g.ip = np.full(g.N_eng, ip_vector[j])
            sub_vect = np.array([alpha_vector[i], beta, p, q, r, da, de, dr])

            V_vect = np.zeros((g.N_eng))

            for k in range(g.N_eng):
                a = np.array([g.xp[k], g.yp[k], g.zp[k]])
                V_vect[k] = (V + np.cross(ang_speed, a)[0])*np.cos((-np.sign(g.yp[k])) * beta + g.wingsweep)*np.cos(alpha_vector[i]+g.alpha_i-g.alpha_0+g.ip[k])

            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
            # CL_tail = AeroForces.Cm_and_CL_tail(V, np.copy(CoefMatrix), x, Tc, atmo, g, PropWing)[1]
            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)




    print('Plot1')
    print(ip_vector)
    print(CL_matrix)

    g.x_offset = original_offset_vector
    """



    #PLOT 2
    #Angles of attack: 4, 8 [°] (two lines)
    #Axe x: Variation of installation angle
    #Axe y: CL,eff,  adding vertical contribution of thrust

    """
    alpha_vector = np.linspace(4, 8, 2)*np.pi/180
    ip_vector = np.linspace(-20, 8, 29)*np.pi/180
    CL_matrix = np.zeros((len(alpha_vector), len(ip_vector)))

    for i in range(len(alpha_vector)):

        for j in range(len(ip_vector)):
            g.ip = np.full(g.N_eng, ip_vector[j])
            sub_vect = np.array([alpha_vector[i], beta, p, q, r, da, de, dr])

            V_vect = np.zeros((g.N_eng))

            for k in range(g.N_eng):
                a = np.array([g.xp[i], g.yp[i], g.zp[i]])
                V_vect[k] = (V + np.cross(ang_speed, a)[0])*np.cos((-np.sign(g.yp[k])) * beta + g.wingsweep)*np.cos(alpha_vector[i]+g.alpha_i-g.alpha_0+g.ip[k])

            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi



            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)

            # Matrix to transform a vector from body reference to aero reference
            Body2Aero_matrix = np.array([[np.cos(alpha_vector[i])*np.cos(beta), np.sin(beta), np.sin(alpha_vector[i])*np.cos(beta)], [-np.cos(alpha_vector[i])*np.sin(beta), np.cos(beta), -np.sin(alpha_vector[i])*np.sin(beta)], [-np.sin(alpha_vector[i]), 0, np.cos(alpha_vector[i])]])

            # Moment and Force of thrust  is obtained in body reference
            Moment = np.zeros((g.N_eng, 3))
            F_thrust_body = np.zeros((g.N_eng, 3))
            for k in range(g.N_eng):
                a = np.array([g.xp[k], g.yp[k], g.zp[k]])
                b = np.array([Fx_vec[k]*np.cos(g.alpha_i - g.alpha_0+g.ip[k]), 0, -Fx_vec[k]*np.sin(g.alpha_i - g.alpha_0+g.ip[k])])
                Moment[k, :] = np.cross(a, b)
                F_thrust_body[k, :] = b
            Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))
            F_thrust_body = np.array((np.sum(F_thrust_body[:, 0]), np.sum(F_thrust_body[:, 1]), np.sum(F_thrust_body[:, 2])))
            Mt = Thrust_moment_body

            # Thrust force is transformed from body to aero reference
            F_thrust_aero = Body2Aero_matrix @ F_thrust_body




            CL_matrix[i, j] = (-F[2] - F_thrust_aero[2] )/ (0.5*rho*V**2 * g.S)

    print('Plot2')
    print(ip_vector)
    print(CL_matrix)

    g.ip = original_ip_vector
    """


    #PLOT 3
    #Angles of attack: 4, 8 [°] (two lines)
    #Axe x: Variation of offset
    #Axe y: CL,wing
    """


    #g.ip = np.full(g.N_eng, -20*np.pi/180)
    alpha_vector = np.linspace(4, 8, 2)*np.pi/180
    offset_vector = np.linspace(0.5*0.254, 6*0.254, 20)
    CL_matrix = np.zeros((len(alpha_vector), len(offset_vector)))

    for i in range(len(alpha_vector)):

        for j in range(len(offset_vector)):
            g.x_offset = np.full(g.N_eng, offset_vector[j])
            sub_vect = np.array([alpha_vector[i], beta, p, q, r, da, de, dr])

            V_vect = np.zeros((g.N_eng))

            for k in range(g.N_eng):
                a = np.array([g.xp[k], g.yp[k], g.zp[k]])
                V_vect[k] = (V + np.cross(ang_speed, a)[0])*np.cos((-np.sign(g.yp[k])) * beta + g.wingsweep)*np.cos(alpha_vector[i]+g.alpha_i-g.alpha_0+g.ip[k])

            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)



    print('Plot3')
    print(offset_vector)
    print(CL_matrix)

    g.x_offset = original_offset_vector
    """




    #PLOT 4 Diameter
    #% Angles of attack: 4, 8 (two lines)
    #% Axe x: Variation of Diameter (with current radius there is no more
    #% space) from D=D to D = 0.25D
    #% Axe y: CL,wing

    """
    alpha_vector = np.linspace(8, 8, 1)*np.pi/180
    nominaldiameter = original_diameter[0]
    diametervector = np.linspace(0.75*nominaldiameter, nominaldiameter, 10)
    CL_matrix = np.zeros((len(alpha_vector), len(diametervector)))

    for i in range(len(alpha_vector)):

        for j in range(len(diametervector)):
            g.Dp = np.full(g.N_eng, diametervector[j])
            sub_vect = np.array([alpha_vector[i], beta, p, q, r, da, de, dr])

            V_vect = np.zeros((g.N_eng))

            for k in range(g.N_eng):
                a = np.array([g.xp[k], g.yp[k], g.zp[k]])
                V_vect[k] = (V + np.cross(ang_speed, a)[0])*np.cos((-np.sign(g.yp[k])) * beta + g.wingsweep)*np.cos(alpha_vector[i]+g.alpha_i-g.alpha_0+g.ip[k])

            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi


            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)



    print('Plot4')
    print(diametervector/nominaldiameter)
    print(CL_matrix)

    g.Dp = original_diameter
    """




    #PLOT 5 Number of engines
    # Varying the number of engines ensuring that all the area is wetted.
    #% Angles of attack: 4, 8 (two lines)
    #% Keep radius
    #% keep total power
    #% Axe x: Variation of number of engines = [18 16 14 12 10 8 6 ]
    #% Axe y: CL,wing

    """
    alpha_vector = np.linspace(8,8, 1)*np.pi/180
    Number_engines = np.array([18, 16, 14, 12, 10, 8, 6])
    CL_matrix = np.zeros((len(alpha_vector), len(Number_engines)))

    for i in range(len(alpha_vector)):

        for j in range(len(Number_engines)):

            g.N_eng = int(Number_engines[j])

            g.yp, RHLP = Distributed_engines(g.N_eng, g.Dp[0])

            g.Dp = np.full(g.N_eng, 2*RHLP)
            g.Sp = g.Dp**2/4*math.pi

            g.xp = np.full(g.N_eng, 10*0.02547)

            g.zp = np.full(g.N_eng, -0.454052)

            g.x_offset = np.full(g.N_eng, 10*0.02547)
            g.ip = np.full(g.N_eng, -7.25/180 * np.pi)

            sub_vect = np.array([alpha_vector[i], beta, p, q, r, da, de, dr])

            V_vect = np.ones((g.N_eng))*V



            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)
   


    print('Plot5')
    print(Number_engines)
    print(CL_matrix)

    """




    #Plot 6 
    #Now for:
    #No flaps
    #alpha = 8
    #dx = 1
    #We change the number of engines.
    #For every number engine we are going to change the diameter between the maximum diameter (for having all the wing blown)
    #and 0.25 times that one. When diameter is the maximum, the position of the engines is fixed.


    """
    alpha = 8*np.pi/180

    Number_engines = np.array([18, 16, 14, 12, 10, 8, 6])

    Diameter = np.linspace(1, 0.25, 7)

    CL_matrix = np.zeros((len(Number_engines), len(Diameter)))

    Is_stall = np.zeros((len(Number_engines), len(Diameter)))
    Is_maxpower = np.zeros((len(Number_engines), len(Diameter)))

    Rtip = 0.5 * 60 * 0.0254  # radius of propeller at tip
    FusR = 0.5*1.266  # max. radius fuselage

    for i in range(len(Number_engines)):

        for j in range(len(Diameter)):

            g.N_eng = int(Number_engines[i])

            g.yp, RHLP = Distributed_engines(g.N_eng, g.Dp[0])


            g.Dp = np.full(g.N_eng, 2*RHLP*Diameter[j])
            g.Sp = g.Dp**2/4*math.pi

            g.xp = np.full(g.N_eng, 10*0.02547)
            g.zp = np.full(g.N_eng, -0.454052)

            g.x_offset = np.full(g.N_eng, 10*0.02547)
            g.ip = np.full(g.N_eng, -7.25/180 * np.pi)



            sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])

            V_vect = np.zeros((g.N_eng))

            for k in range(g.N_eng):
                a = np.array([g.xp[k], g.yp[k], g.zp[k]])
                V_vect[k] = (V + np.cross(ang_speed, a)[0])*np.cos((-np.sign(g.yp[k])) * beta + g.wingsweep)*np.cos(alpha+g.alpha_i-g.alpha_0+g.ip[k])

            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)



    print('Plot6')
    print(Number_engines)
    print(Diameter)
    print(CL_matrix)
    """




    #Plot 7
    #Now for:
    #No flaps
    #alpha = 8
    #dx = 1
    #We change the number of engines
    #For every number engine we are going to change the diameter between the maximum diameter (for having all the wing blown)
    #and 0.6 times that one. When diameter is the maximum, the position of the engines is fixed. If we reduce the diameter
    #Then we will approach the engines to the one close to the fuselage.



    alpha = 8*np.pi/180

    Number_engines = np.array([18, 16, 14, 12, 10, 8, 6])

    Diameter = np.linspace(1, 0.8, 3)

    CL_matrix = np.zeros((len(Number_engines), len(Diameter)))

    Is_stall = np.zeros((len(Number_engines), len(Diameter)))
    Is_maxpower = np.zeros((len(Number_engines), len(Diameter)))

    Rtip = 0.5 * 60 * 0.0254  # radius of propeller at tip
    FusR = 0.60198  # max. radius fuselage

    for i in range(len(Number_engines)):

        for j in range(len(Diameter)):

            g.N_eng = int(Number_engines[i])

            g.yp, RHLP = Distributed_engines(g.N_eng, g.Dp[0])


            g.Dp = np.full(g.N_eng, 2*RHLP*Diameter[j])
            g.Sp = g.Dp**2/4*math.pi

            g.xp = np.full(g.N_eng, 10*0.02547)
            g.zp = np.full(g.N_eng, -0.454052)

            g.x_offset = np.full(g.N_eng, 10*0.02547)
            g.ip = np.full(g.N_eng, -7.25/180 * np.pi)

            # We need to redifine yp
            if Diameter[j] == 1:
                pass
            else:
                min_y = (FusR + RHLP*Diameter[j])
                max_y = (min_y + (0.5*g.N_eng-1)*2*RHLP*Diameter[j])

                yp = np.linspace(min_y, max_y, int(g.N_eng/2))
                g.yp = np.hstack((np.flip(-yp), yp))


            sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])

            V_vect = np.ones((g.N_eng))*V

            Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
            Fx = np.sum(Fx_vec)

            # Convert thrust in Tc for patterson
            Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

            F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
            CL_matrix[i, j] = -F[2] / (0.5*rho*V**2 * g.S)


            x = np.array([alpha, 0, 0, 0, 0, 0, 0, 0, 0])
            fixtest = np.array([V, 0, 0, 0])
            for h in range(int(g.N_eng)):
               x = np.append(x, dx)
            #printx(x, fixtest, atmo, g, PropWing)



    print('Plot7')
    print(Number_engines)
    print(Diameter)
    print(CL_matrix)





    #Plot 8
    #The same as before, but now we vary alpha.


    """
    alpha_vector = np.linspace(-5,20,26)*np.pi/180

    Number_engines = np.array([18, 16, 14, 12, 10, 8, 6])

    Diameter = np.linspace(1, 1, 1)

    CL_matrix = np.zeros((len(Number_engines), len(Diameter), len(alpha_vector)))


    Rtip = 0.5 * 60 * 0.0254  # radius of propeller at tip
    FusR = 0.60198  # max. radius fuselage

    print('Plot8')
    for i in range(len(Number_engines)):

        print(Number_engines[i])

        for j in range(len(Diameter)):
            print(Diameter[j])

            for counter in range(len(alpha_vector)):


                  alpha = alpha_vector[counter]
                  print(alpha*180/np.pi)

                  g.N_eng = int(Number_engines[i])

                  g.yp, RHLP = Distributed_engines(g.N_eng, g.Dp[0])


                  g.Dp = np.full(g.N_eng, 2*RHLP*Diameter[j])
                  g.Sp = g.Dp**2/4*math.pi

                  g.xp = np.full(g.N_eng, 10*0.02547)
                  g.zp = np.full(g.N_eng, -0.454052)

                  g.x_offset = np.full(g.N_eng, 10*0.02547)
                  g.ip = np.full(g.N_eng, -7.25/180 * np.pi)

                  # We need to redifine yp
                  if Diameter[j] == 1:
                      pass
                  else:
                      min_y = (FusR + RHLP*Diameter[j])
                      max_y = (min_y + (0.5*g.N_eng-1)*2*RHLP*Diameter[j])

                      yp = np.linspace(min_y, max_y, int(g.N_eng/2))
                      g.yp = np.hstack((np.flip(-yp), yp))


                  sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])

                  V_vect = np.ones((g.N_eng))*V

                  Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
                  Fx = np.sum(Fx_vec)

                  # Convert thrust in Tc for patterson
                  Tc = Fx_vec/(2*rho*g.Sp*V**2)  # For turning dimensionless use V, V_vect has already been used for calculating FXi

                  F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
                  CL_matrix[i, j, counter] = -F[2] / (0.5*rho*V**2 * g.S)
                  print(CL_matrix[i,j,counter])







    print('FINAL PRINT')
    print(Number_engines)
    print(Diameter)
    print(alpha_vector*180/np.pi)
    print(CL_matrix)
    """




    print('End')


    return CL_matrix







