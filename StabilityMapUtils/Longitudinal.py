import numpy as np
import math
from scipy.optimize import minimize
from scipy.optimize import least_squares
import sys
from StabilityMapUtils import AeroForces




def Long_Equilibrium(Coef_base, atmospher, g, PW, vars) :



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
          objective function is e.V_min) (OPTIMIZATION PERFORMED).
        * By defining V and dx    the function will find the set of (alpha de gamma theta) to accomplish equilibrium.
        * By defining V and gamma the function will find the set of (alpha de dx    theta)  to accomplish equilibrium.
    """

    MaxIter = 1000
    tolerance = 1e-10

    alpha = vars[0]
    gamma = vars[1]
    V = vars[2]
    de = vars[3]
    dx = vars[4]
    theta = vars[5]




    Variables = {'alpha': alpha, 'theta': theta, 'de': de, 'dx': dx , 'V': V, 'gamma': gamma }
    Bounds = {'alpha': (-5*math.pi/180, 20*math.pi/180), 'theta': (-30/180*math.pi, 30/180*math.pi), 'de': (-23/180*math.pi, 13/180*math.pi), 'dx': (1e-9, 1), 'V': (0, 70), 'gamma': (-30/180*math.pi, 30/180*math.pi)}
    Initial_guess = {'alpha': 0.10615368859616711, 'theta': 0.10615368859618061, 'de': -0.0924842169820762, 'dx': 0.3122270558364951, 'V': 70, 'gamma': 0}

    x0 = []
    fixtest = []
    bnds1 = []
    bnds2 = []

    xorder = []
    fixtestorder = []

    for key, value in Variables.items():

        if type(value) != str:
               fixtest.append(value)
               fixtestorder.append(key)

        else:
                x0.append(Initial_guess[key])
                xorder.append(key)
                bnds1.append(Bounds[key][0])
                bnds2.append(Bounds[key][1])



    fixtest = np.array(fixtest)
    x0 = np.array(x0)
    bnds1 = np.array(bnds1)
    bnds2 = np.array(bnds2)

    diccons = (np.copy(fixtest), xorder, fixtestorder, np.copy(Coef_base), atmospher, g, PW)




    if len(x0) == 4:   # Solving acotated problem
        bnds = (bnds1, bnds2)
        # V = 70 , gamma = 0
        # xo = alpha = 0.10615368859616711 theta = 0.10615368859618061 , dx = 0.3122270558364951 , de = -0.0924842169820762
        k = least_squares(Order, x0, args=diccons, bounds=bnds, max_nfev = 2000, ftol=1e-8, xtol=1e-8)

    elif len(x0) >= 5:  # Solving optimization problem. An objective function shall be given.

        bnds = []

        for i in range(len(x0)):
            bnds.append((bnds1[i], bnds2[i]))

        bnds = tuple(bnds)




        k = minimize(V_min, x0, args=diccons, bounds=bnds,
                     constraints={'type': 'eq', 'fun': Order, 'args': diccons},
                     options={'maxiter': MaxIter, 'disp': True}, tol=tolerance)

    else:
        print("Problem may not be well defined")
        sys.exit()

    # check if constraints are validated
    constraints_calc = Order(k.x, *diccons)
    print("\nConstraints")





    i = -1
    for key in xorder:
        i = i+1
        if key == 'V':
            V = k.x[i]
        elif key == 'alpha':
            alpha = k.x[i]
        elif key == 'de':
            de = k.x[i]
        elif key == 'dx':
            dx = k.x[i]
        elif key == 'theta':
            theta = k.x[i]
        elif key == 'gamma':
             gamma = k.x[i]
    i = -1
    for key in fixtestorder:
        i = i+1
        if key == 'V':
            V = fixtest[i]
        elif key == 'alpha':
            alpha = fixtest[i]
        elif key == 'de':
            de = fixtest[i]
        elif key == 'dx':
            dx = fixtest[i]
        elif key == 'theta':
            theta = fixtest[i]
        elif key == 'gamma':
             gamma = fixtest[i]

    x = np.concatenate((np.array([alpha, 0, 0, 0, 0, theta, 0, de, 0]), np.full(g.N_eng, dx)))
    fixtest = np.array([V, 0, gamma, 0])

    # for coherance shall be
    #       x = [alpha, p, q ,r, phi, theta, da, de, dr, dx]
    # fixtest = [V, beta, gamma, omega]

    return k, x, fixtest,  (fixtest, np.copy(Coef_base), atmospher, g, PW)






def Order(a, b, xorder, fixtestorder,CoefMatrix, atmo, g, PropWing):

    """
    You need at the end a vector like this
    -x =[V, alpha, theta, de, dx, gamma]


    """

    x = np.zeros(6)
    i = -1

    for key in xorder:
        i = i+1

        if key == 'V':
            x[0] = a[i]

        elif key == 'alpha':
            x[1] = a[i]

        elif key == 'theta':
            x[2] = a[i]

        elif key == 'de':
            x[3] = a[i]

        elif key == 'dx':
            x[4] = a[i]

        elif key == 'gamma':
            x[5] = a[i]


    i = -1
    for key in fixtestorder:
        i = i+1

        if key == 'V':
            x[0] = b[i]

        elif key == 'alpha':
            x[1] = b[i]

        elif key == 'theta':
            x[2] = b[i]

        elif key == 'de':
            x[3] = b[i]

        elif key == 'dx':
            x[4] = b[i]

        elif key == 'gamma':
            x[5] = b[i]


    A = Long_equations(x, CoefMatrix, atmo, g, PropWing)

    return A












def V_min(a, b, xorder, fixtestorder,CoefMatrix, atmo, g, PropWing):


    rho = atmo[1]

    # --- Now prepare variables for equations ---
    i = -1
    for key in xorder:
        i = i+1
        if key == 'V':
            V = a[i]
        elif key == 'alpha':
            alpha = a[i]
        elif key == 'de':
            de = a[i]
        elif key == 'dx':
            dx = a[i]
    i = -1
    for key in fixtestorder:
        i = i+1
        if key == 'V':
            V = b[i]
        elif key == 'alpha':
            alpha = b[i]
        elif key == 'de':
            de = b[i]
        elif key == 'dx':
            dx = b[i]

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
    Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)
    Fx = np.sum(Fx_vec)

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)

    f = ((9.81*g.m - Fx * np.sin(g.alpha_i + g.alpha_0+g.ip + alpha))/(np.abs(F[2])/V**2))**0.5

    return f








def Long_equations(x, CoefMatrix, atmo, g, PropWing):
    """function defining constraints for speed minimization in longitudinal
    inputs:
        -x =[V, alpha, theta, delta_e, delta_i]
        length of x except the propulsion levels is 8
        -fix = [gamma, beta, p, q, r, phi, da, dr]

        gamma = beta = p = q = r = phi = da = dr = 0 as we are in LONGITUDINAL equilibrium

        4 equations (2 forces, 1 moment, theta = alpha + gamma)
    """


    rho = atmo[1]

    # --- Now prepare variables for equations ---
    V = x[0]
    alpha = x[1]
    theta = x[2]
    de = x[3]
    dx = x[4]
    gamma = x[5]


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



    Fx_vec = g.Thrust(np.full(g.N_eng, dx), V_vect, atmo)
    Fx = np.sum(Fx_vec)



    #Matrix to transform a vector from body reference to aero reference
    Body2Aero_matrix = np.array([[np.cos(alpha)*np.cos(beta), np.sin(beta), np.sin(alpha)*np.cos(beta)], [-np.cos(alpha)*np.sin(beta), np.cos(beta), -np.sin(alpha)*np.sin(beta)], [-np.sin(alpha), 0, np.cos(alpha)]])

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






    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)                                                                                       #For adimension V, has already been used for calculating FXi

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)


    #F gives out aerodinamical forces in aero axis: Drag, lateral force and lift and moments
    # Does not give out X,Y,Z
    Aero2Body_matrix = np.transpose(Body2Aero_matrix)

    F_aero_bodyref = Aero2Body_matrix @ F[0:3]










    A = np.zeros(4)

    """
    Version 1 
    Equations in aero frame. 
    By using this formulation you get directly V_dot, alpha_dot, gamma_dot, but the equations
    expressed like this (divided by g.m or g.Iy) are way smaller and its more easier to accomplish with 
    A = 0 with same tolerance. In other words, if you write them like in version 2 you are decrasing tolerance.
    
    A[0] = +(F[0] + F_thrust_aero[0])/g.m
    A[1] = 9.81/V + (F[2] + F_thrust_aero[2])/(g.m*V)
    A[2] = (Mt[1] + F[4])/g.Iy
    A[3] = alpha + gamma - theta
    """

    """
    Version 2. 
    Equations in aero frame. Keep in mind there are not V_dot, alpha_dot or gamma_dot
    A[0] = +(F[0] + F_thrust_aero[0]) - 9.81*g.m*np.sin(gamma)
    A[1] = 9.81*g.m*np.cos(gamma) + (F[2] + F_thrust_aero[2])
    A[2] = (Mt[1] + F[4])
    A[3] = alpha + gamma - theta
    """





    A[0] = +(F_aero_bodyref[0] + F_thrust_body[0]) - 9.81*g.m*np.sin(theta)
    A[1] = 9.81*g.m*np.cos(theta) + (F_aero_bodyref[2] + F_thrust_body[2])
    A[2] = (Mt[1] + F[4])
    A[3] = alpha + gamma - theta



    return A
