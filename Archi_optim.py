"""
File for performing the optimization of the architecture of the propellers for minimim V_stall
We tried to minimize the stall speed while changing the architecture. In construction...


@david.planas


[ 0.28837467 30.         -0.81914941  0.03748512  0.08398621  0.10067369
  0.10961772  0.11126013  0.10916654  0.10299067  0.08656456  0.3454908
  0.3454908   0.3454908   0.3454908   0.3454908   0.3454908   0.3454908
  0.3454908 ] DERNIERE ITERATION;
"""


import numpy as np
import math
from scipy.optimize import minimize
from scipy.optimize import least_squares
import sys
from StabilityMapUtils import AeroForces



def Archi_optim(Coef_base, atmospher, g, PW):

     MaxIter = 1000
     tolerance = 1e-10
     g.N_eng = 16      # N_eng : 16 14 12 10 8 6

     b = 9.642  # wingspan
     Rtip = 0.5 * 60 * 0.0254  # radius of propeller at tip
     FusR = 0.60198  # max. radius fuselage

     Dporiginal = g.Dp[0]
     Dpbounds = (0.6*Dporiginal, 2*Dporiginal),
     dx = (1e-9, 1),
     dx0 = 0.8

     x0 = [0, 10, 0]
     bounds = ((-5*np.pi/180, 20*np.pi/180), (10, 30), (-100*np.pi/180, 100*np.pi/180))

     for i in range(int(g.N_eng/2)):
         bounds = bounds + dx
         x0.append(dx0)

     for j in range(int(g.N_eng/2)):
         bounds = bounds + Dpbounds
         x0.append((0.5*b - Rtip - FusR)/(g.N_eng/2))

     x0 = np.array(x0)
     diccons = (np.copy(Coef_base), atmospher, g, PW)

     cons = ({'type': 'eq', 'fun': Constraints, 'args': diccons},
             {'type': 'ineq', 'fun': lambda x: (0.5*b - Rtip) - (sum(x[3:])+FusR)})

     k = minimize(V_min,  # Function to minimize
                  x0,  # Initial value of the vector of variables to vary. Numpy array.
                  args=diccons,  # Extra arguments for the function to minimize.
                  bounds=bounds,  # Possible bounds for the vector of variables to vary. Tuple
                  constraints=cons,  # Defining the constraints. Type: eq for equality. fun: the function with constraints. args: the arguments to the function.
                  options={'maxiter': MaxIter, 'disp': True}, tol=tolerance)  # Options for the optimizer.

     return k











def V_min(x,CoefMatrix, atmo, g, PropWing):


     FusR = 0.60198

     g.Dp = np.hstack((x[int(-g.N_eng/2):], np.flip(x[int(-g.N_eng/2):])))
     g.Sp = g.Dp**2/4*math.pi

     g.xp = np.full(g.N_eng, 10*0.02547)
     g.zp = np.full(g.N_eng, -0.454052)

     g.x_offset = np.full(g.N_eng, 10*0.02547)
     g.ip = np.full(g.N_eng, -7.25/180 * np.pi)

     yp = np.zeros(int(0.5*g.N_eng))

     for i in range(len(yp)):
         if i == 0:
             yp[-1] = FusR + g.Dp[-1]/2
         else:
             yp[-1-i] = yp[-i] + g.Dp[-i]/2 + g.Dp[-i-1]/2

     g.yp = np.hstack((np.flip(-yp), yp))

     alpha = x[0]
     V = x[1]
     de = x[2]
     dx = np.hstack((x[-g.N_eng:-int(g.N_eng/2)], np.flip(x[-g.N_eng:-int(g.N_eng/2)])))

     rho = atmo[1]
     p = 0
     q = 0
     r = 0
     beta = 0
     da = 0
     dr = 0



     # --- Compute aerodynamic forces ---
     #here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
     sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])  # rudder is allowed

     #Thrust forces and moments
     V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp
     Fx_vec = g.Thrust(dx, V_vect, atmo)
     # convert thrust in Tc for patterson
     Tc = Fx_vec/(2*rho*g.Sp*V**2)
     Fx = np.sum(Fx_vec)

     F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)




     # Moment and Force of thrust  is obtained in body reference
     F_thrust_body = np.zeros((g.N_eng, 3))
     for i in range(g.N_eng):
         F_thrust_body[i, :] = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip[i]), 0, -Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip[i])])
     F_thrust_body = np.array((np.sum(F_thrust_body[:, 0]), np.sum(F_thrust_body[:, 1]), np.sum(F_thrust_body[:, 2])))


     f = ((9.81*g.m - F_thrust_body[0]*np.sin(alpha) - F_thrust_body[2]*np.cos(alpha))/(np.abs(F[2])/V**2))**0.5

     # f = ((9.81*g.m - Fx * np.sin(g.alpha_i + g.alpha_0+g.ip + alpha))/(np.abs(F[2])/V**2))**0.5 i did ICAS with this, I think sign of alpha_0 is wrong

     return f










def Constraints(x, CoefMatrix, atmo, g, PropWing):
    """function defining constraints for speed minimization in longitudinal
    inputs:
        -x =[V, alpha, theta, delta_e, delta_i]
        length of x except the propulsion levels is 8
        -fix = [gamma, beta, p, q, r, phi, da, dr]

        gamma = beta = p = q = r = phi = da = dr = 0 as we are in LONGITUDINAL equilibrium

        4 equations (2 forces, 1 moment, theta = alpha + gamma)
    """

    print(x)




    FusR = 0.60198

    g.Dp = np.hstack((x[int(-g.N_eng/2):], np.flip(x[int(-g.N_eng/2):])))
    g.Sp = g.Dp**2/4*math.pi

    g.xp = np.full(g.N_eng, 10*0.02547)
    g.zp = np.full(g.N_eng, -0.454052)

    g.x_offset = np.full(g.N_eng, 10*0.02547)
    g.ip = np.full(g.N_eng, -7.25/180 * np.pi)

    yp = np.zeros(int(0.5*g.N_eng))

    for i in range(len(yp)):
        if i == 0:
            yp[-1] = FusR + g.Dp[-1]/2
        else:
            yp[-1-i] = yp[-i] + g.Dp[-i]/2 + g.Dp[-i-1]/2

    g.yp = np.hstack((np.flip(-yp), yp))



    rho = atmo[1]

    # --- Now prepare variables for equations ---
    alpha = x[0]
    V = x[1]
    theta = x[0]
    de = x[2]
    dx = np.hstack((x[-g.N_eng:-int(g.N_eng/2)], np.flip(x[-g.N_eng:-int(g.N_eng/2)])))
    gamma = 0

    p = 0
    q = 0
    r = 0
    beta = 0
    da = 0
    dr = 0


    # --- Compute aerodynamic forces ---
    #here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect = np.array([alpha, beta, p, q, r, da, de, dr])  # rudder is allowed


    #Thrust forces and moments

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.yp)) * beta + g.wingsweep) - r * g.yp


    Fx_vec = g.Thrust(dx, V_vect, atmo)
    Fx = np.sum(Fx_vec)

    #Matrix to transform a vector from body reference to aero reference
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

    # Moment and Force of thrust  is obtained in body reference
    Moment = np.zeros((g.N_eng, 3))
    F_thrust_body1 = np.zeros((g.N_eng, 3))
    for i in range(g.N_eng):
        a = np.array([g.xp[i], g.yp[i], g.zp[i]])
        b = np.array([Fx_vec[i]*np.cos(g.alpha_i - g.alpha_0+g.ip[i]), 0, -Fx_vec[i]*np.sin(g.alpha_i - g.alpha_0+g.ip[i])])
        Moment[i, :] = np.cross(a, b)
        F_thrust_body1[i, :] = b
    Thrust_moment_body = np.array((np.sum(Moment[:, 0]), np.sum(Moment[:, 1]), np.sum(Moment[:, 2])))
    F_thrust_body = np.array((np.sum(F_thrust_body1[:, 0]), np.sum(F_thrust_body1[:, 1]), np.sum(F_thrust_body1[:, 2])))

    Mt = Thrust_moment_body

    # Thrust force is transformed from body to aero reference
    F_thrust_aero = Body2Aero_matrix @ F_thrust_body

    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)                                                                                       #For adimension V, has already been used for calculating FXi

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)

    #F gives out aerodinamical forces in aero axis: Drag, lateral force and lift and moments
    # Does not give out X,Y,Z
    Aero2Body_matrix = np.transpose(Body2Aero_matrix)

    F_aero_bodyref = Aero2Body_matrix @ F[0:3]

    A = np.zeros(3)

    A[0] = +(F_aero_bodyref[0] + F_thrust_body[0]) - 9.81*g.m*np.sin(theta)
    A[1] = 9.81*g.m*np.cos(theta) + (F_aero_bodyref[2] + F_thrust_body[2])
    A[2] = (Mt[1] + F[4])


    return A

