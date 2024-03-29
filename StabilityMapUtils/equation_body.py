"""
12/05/2022

Module defining non linear flight equations. An attempt to write them in body reference system.

@author: e.nguyen-van
david.planas-andres
"""
import numpy
import numpy as np
import math
from numpy.linalg import inv
from StabilityMapUtils import AeroForces





def Constraints_DEP(x, fix, CoefMatrix, atmo, g, PropWing):
    """function defining constraints for power minimization
    inputs:
    -x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
    x is the state to determine
    length of x except the propulsion levels is 8
    -fix = [V, beta, gamma, omega]
    fix is the vector of parameters whom are fixed by the user

    """


    rho = atmo[1]

    # --- Now prepare variables for equations ---
    V = fix[0]
    alpha = x[0]
    beta = fix[1]
    gamma = fix[2]
    omega = fix[-1]
    p = x[1]
    q = x[2]
    r = x[3]
    phi = x[4]
    theta = x[5]
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
    F[0:3] = np.transpose(Body2Aero_matrix) @ F[0:2]
    F[3:6] = np.transpose(Body2Aero_matrix) @ F[0:2]


    # Transformation of aerodynamic speed from aero to body reference
    [u,v,w] = np.transpose(Body2Aero_matrix) @ np.concatenate(([V], [0], [0]))





    A = np.zeros(10+g.inop)


    A[0] = (1/g.m)*(F_thrust_body[0] + F[0]) - 9.81*np.sin(theta) +r*v - q*w
    A[1] = (1/g.m)*(F_thrust_body[0] + F[1]) + 9.81*np.cos(theta)*np.sin(phi) - r*u + p*w
    A[2] = (1/g.m)*(F_thrust_body[2] + F[2]) + 9.81*np.cos(theta)*np.cos(phi) + q*u - p*v
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















def fobjectivePropWingInterac(x, fix, rho, g):

    Dx = x[-g.N_eng:]
    MeanDx = np.mean(Dx)
    stdDx = np.std(Dx)

    if g.hangar['aircraft'] == 'ATR72':
         #    Power=np.sum(x[-g.N_eng:])*2*g.P_var/float(g.N_eng)*rho/1.225/1000000
         return MeanDx+stdDx

    elif g.hangar['aircraft'] == 'DECOL':
         return MeanDx*0.5+stdDx*0.5















def Jac_DEP(x, fix, CoefMatrix, atmo, g, PropWing, h):
    # function to compute the jacobian at a steady state
    # the function is hard coded inside
    # inputs :
    #       -x : steady state vector
    #       -fixtuple : tuple of (fixed param, function param)
    #       -h : step to compute derivative

    nfx=9 # number of equations for flight analysis (V, beta, alpha, p, q, r, phi, theta, gamma)
    # As gamma is a parameter in the flight equation (gamma_dot not computed),
    # the vector of accelerations is : [V,beta,alpha,p,q,r,phi,theta] = nfx-1

    step_vec=x*h

    for i in range(len(step_vec)):
        # check for zeros
        if step_vec[i]<1e-4:
            step_vec[i]=0.001

#    fx=Constraints_DEP(x, *fixtuple)

    dx=np.zeros((nfx-1,len(x)+3))
    fixtuple=(fix, CoefMatrix, atmo, g, PropWing)

    # compute derivative using centered difference
    #Accelerations due to a small change in velocity
    fix_plus=fix+np.append([fix[0]*h/2.0],np.zeros((len(fix)-1)))
    fix_minus=fix-np.append([fix[0]*h/2.0],np.zeros((len(fix)-1)))

    tuple_plus=(fix_plus, CoefMatrix, atmo, g, PropWing)
    tuple_minus=(fix_minus, CoefMatrix, atmo, g, PropWing)

    diff=(Constraints_DEP(x,*tuple_plus)-Constraints_DEP(x,*tuple_minus))/(fix[0]*h)
    dx[:,0]=diff[0:nfx-1]

    #Accelerations due to a small change in side-slip
    beta_step=np.zeros((len(fix)))
    beta_step[1]=h/2
    fix_plus=fix+beta_step
    fix_minus=fix-beta_step

    tuple_plus=(fix_plus, CoefMatrix, atmo, g, PropWing)
    tuple_minus=(fix_minus, CoefMatrix, atmo, g, PropWing)

    diff=(Constraints_DEP(x,*tuple_plus)-Constraints_DEP(x,*tuple_minus))/(beta_step[1]*2)
    dx[:,1]=diff[0:nfx-1]

    #Accelerations due to a small change in gamma
    gamma_step=np.zeros((len(fix)))
    gamma_step[2]=h/2
    fix_plus=fix+gamma_step
    fix_minus=fix-gamma_step

    tuple_plus=(fix_plus, CoefMatrix, atmo, g, PropWing)
    tuple_minus=(fix_minus, CoefMatrix, atmo, g, PropWing)

    diff=(Constraints_DEP(x,*tuple_plus)-Constraints_DEP(x,*tuple_minus))/(gamma_step[2]*2)
    dx[:,2]=diff[0:nfx-1]

    #now all acceleration due to a small of each variables in x
    for j in range(len(x)):
        activex=np.zeros((len(x)))
        activex[j]=1
        dfx=(Constraints_DEP(x+activex*step_vec/2,*fixtuple)-Constraints_DEP(x-activex*step_vec/2,*fixtuple))/np.dot(activex,step_vec)
        dx[:,j+3]=dfx[0:nfx-1]

    # optionally decouple matrix
    return dx