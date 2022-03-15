# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 17:01:52 2017

Module defining non linear flight equations

@author: e.nguyen-van
         david.planas-andres
"""

import numpy as np
import math
from numpy.linalg import inv
from StabilityMapUtils import AeroForces
# functions for standard earth equation

def X(x, CoefMatrix,Velocities,rho,g):
    # x : state [u,v,w,p,q,r,phi,theta,dfl,da,de,dr,dx]
    # arg : tuple of additional elements, respect order below:
#    CoefMatrix=arg[0]
#    Velocities = arg[1]
#    rho = arg[2]
#    g = arg[3]


    con=np.array([80,0/180*math.pi,0,0,0,0,0,0,0,0,0,0])
    if x[7]<-30/180*math.pi:
        return np.ones(12)*1000-np.ones(12)*x[7]*1000
    elif x[7]>30/180*math.pi:
        return np.ones(12)*1000+np.ones(12)*x[7]*1000
    
    if x[-1]<0:
        return np.ones(12)*1000-np.ones(12)*x[-1]*10000
    elif x[-1]>1:
        return np.ones(12)*1000+np.ones(12)*x[-1]*1000
    
    # need to do some operation before computing the forces
    V=math.sqrt(np.dot(x[0:3],x[0:3]))
    sub_vect=np.array([math.atan(x[2]/x[0]), math.asin(x[1]/V)])
    sub_vect=np.append(sub_vect,x[3:6])
    #(alpha, beta, p, q, r, da, dr, de, dx)
    sub_vect=np.append(sub_vect,[x[-4],x[-2],x[-3],x[-1]])

    
    F=AeroForces.CalcForces(V, sub_vect, CoefMatrix, Velocities, rho)
#    print("Printing forces :")
#    print(F)
    #all input taken into account in force computation
    
    vel_vec=x[0:3]
    p=x[3]
    q=x[4]
    r=x[5]
    xdot=-np.cross(np.array([p,q,r]),vel_vec) + 9.81*np.array([-math.sin(x[7]), math.cos(x[7])*math.sin(x[6]), math.cos(x[7])*math.cos(x[6])])+F[0:3]/g.m
#    print(np.cross(np.array([p,q,r]),vel_vec))
    I_inv=np.array([[g.Iz/g.Ix, 0, g.Ixz/g.Ix],
    [0, 1/g.Iy*(g.Iz-g.Ixz**2/g.Ix), 0],
    [g.Ixz/g.Ix, 0, 1]])/(g.Iz-g.Ixz**2/g.Ix)
   
    rot=np.array([p,q,r])
    I=np.array([[g.Ix, 0, -g.Ixz],[0,g.Iy,0],[-g.Ixz,0,g.Iz]]) 
    Mvect=F[3:6]-np.cross(rot,np.array([g.hp,0,0]))-np.cross(rot, np.dot(I,rot))
    
    Mdot=np.dot(I_inv,Mvect);
    xdot=np.append(xdot,Mdot)
    
    xadd=np.empty((0,6))
    xadd=V-con[0]
    xadd=np.append(xadd, [math.asin(x[1]/V)-con[1]])
    xadd=np.append(xadd, [x[3:6]-con[2:5]])
    xadd=np.append(xadd, [x[6]-con[5]])
#    print("xadd:")
#    print(xadd)
#    print("xdot:")
#    print(xdot)
#    
    xreturn=np.append([xdot],[xadd])
    
    return xreturn



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
    V=fix[0]
    alpha=x[0]
    beta=fix[1]
    gamma=fix[2]
    omega=fix[-1]
    p=x[1]
    q=x[2]
    r=x[3]
    phi=x[4]
    theta=x[5]
    I=np.array([ [g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz] ])
    
    # --- Compute aerodynamic forces ---
    #here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect=np.array([alpha,beta,p,q,r])
    if g.nofin==False:
        sub_vect=np.append(sub_vect,[x[6],x[7],x[8]]) # rudder is allowed
    else:
        sub_vect=np.append(sub_vect,[x[6],x[7]]) # no fin allowed, default case



    #Thrust forces and moments

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng


    Fx_vec=g.Thrust(x[-g.N_eng:],V_vect)
    Fx = np.sum(Fx_vec)




    Moment= np.zeros((g.N_eng,3))
    for i in range(g.N_eng):
        a= np.array([ g.x_m , g.PosiEng[i] , g.z_m])
        b=np.array([ Fx_vec[i]*np.cos(g.alpha_i + g.alpha_0+g.ip)  ,  0  ,  -Fx_vec[i]*np.sin(g.alpha_i + g.alpha_0+g.ip)  ])
        Moment[i,:] = np.cross(a,b)
    Thrust_moment_body_axis =np.array(( np.sum(Moment[:,0]), np.sum(Moment[:,1]) , np.sum(Moment[:,2]) ) )

    Body2Aero_matrix = np.array([   [np.cos(alpha)*np.cos(beta), np.sin(beta) , np.sin(alpha)*np.cos(beta) ], [ -np.cos(alpha)*np.sin(beta) , np.cos(beta) , -np.sin(beta)*np.sin(beta) ] , [ -np.sin(alpha), 0   , np.cos(alpha)  ]])

    Trust_moment_aero_axis =    Body2Aero_matrix @  Thrust_moment_body_axis

    Mt = Trust_moment_aero_axis





    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)                                                                                       #For adimension V, has already been used for calculating FXi
    
    F=AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)


    #F gives out aerodinamical forces in aero axis: Drag, lateral force and lift and moments
    # Does not give out X,Y,Z





#     Now sum up the constraints:
    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi) 
    
    A=np.zeros(10+g.inop)
    """
    A created with 10 + number of inoperative engines
    
    A0 = x
    A1 = y
    A2 = z
    A3 = l
    A4 = m
    A5 = n
    A6 = phi
    A7 = theta
    A8 = gamma
    A9 = Omega
    """
    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx*np.cos(alpha+g.alpha_i+g.alpha_0+g.ip)*np.cos(beta)/g.m
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m*V)
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))-Fx*np.sin(alpha+g.alpha_i+g.alpha_0+g.ip)/(g.m*V*np.cos(beta))
    A[3:6]=np.dot(inv(I), np.array([Mt[0],Mt[1],Mt[2]])+F[3:6]-np.cross(np.array([p,q,r]),np.dot(I,np.array([p,q,r]))))
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7]=q*math.cos(phi) -r*math.sin(phi)
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)


    
    for i in range(g.inop):
        A[-1-i]=x[-1-i]                                                                                                 #The inoperative engine are the last ones (right wing). Its value is minimized (to zero)
    
    if g.hangar['version']=='original':                                                                                 #For obligating all the engines to have the same thrust
        #no DEP with original twin or N engines; all engines have the same thrust
        D=np.copy(A)
        for i in range(g.N_eng-g.inop-1):
            AAd=x[-g.N_eng]-x[-g.N_eng+i+1]
            D=np.append(D,[AAd])
        return D
    else:
        return A












def Constraints_minimum_alpha(x, fix, CoefMatrix, atmo, g, PropWing):
    """function defining constraints for alpha minimization
        inputs:
            -x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r,gamma,V]
            x is the state to determine
            length of x is 10
            -fix = [beta, omega, delta_i]
            fix is the vector of parameters whom are fixed by the user, beta (sideslip),
            omega (turn parameter) and  delta_i  (the position of the throttle in percentage of each engine)

    """

    rho = atmo[1]

    V=x[-1]
    alpha=x[0]
    beta=fix[1]
    gamma=x[-2]
    omega=fix[2]
    p=x[1]
    q=x[2]
    r=x[3]
    phi=x[4]
    theta=x[5]
    I=np.array([ [g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz] ])


    # --- Compute aerodynamic forces ---
    #here subvector  must be : (alpha, beta, p, q, r, da, de,dr)
    sub_vect=np.array([alpha,beta,p,q,r])
    if g.nofin==False:
        sub_vect=np.append(sub_vect,[x[6],x[7],x[8]]) # rudder is allowed
    else:
        sub_vect=np.append(sub_vect,[x[6],x[7]]) # no fin allowed, default case


    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng

    Fx_vec=g.Thrust(fix[-g.N_eng:],V_vect)
    Fx = np.sum(Fx_vec)



    Moment= np.zeros((g.N_eng,3))
    for i in range(g.N_eng):
        a= np.array([ g.x_m , g.PosiEng[i] , g.z_m])
        b=np.array([ Fx_vec[i]*np.cos(g.alpha_i + g.alpha_0+g.ip)  ,  0  ,  -Fx_vec[i]*np.sin(g.alpha_i + g.alpha_0+g.ip)  ])
        Moment[i,:] = np.cross(a,b)
    Thrust_moment_body_axis =np.array(( np.sum(Moment[:,0]), np.sum(Moment[:,1]) , np.sum(Moment[:,2]) ) )

    Body2Aero_matrix = np.array([   [np.cos(alpha)*np.cos(beta), np.sin(beta) , np.sin(alpha)*np.cos(beta) ], [ -np.cos(alpha)*np.sin(beta) , np.cos(beta) , -np.sin(beta)*np.sin(beta) ] , [ -np.sin(alpha), 0   , np.cos(alpha)  ]])

    Trust_moment_aero_axis =    Body2Aero_matrix @  Thrust_moment_body_axis

    Mt = Trust_moment_aero_axis




    # convert thrust in Tc for patterson
    Tc = Fx_vec/(2*rho*g.Sp*V**2)

    F=AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)








    #     Now sum up the constraints:
    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi)

    A=np.zeros(10+g.inop)
    """
    A0 = x
    A1 = y
    A2 = z
    A3 = l
    A4 = m
    A5 = n
    A6 = phi
    A7 = theta
    A8 = gamma
    A9 = Omega
    """
    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx*np.cos(alpha+g.alpha_i+g.alpha_0+g.ip)*np.cos(beta)/g.m
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m*V)
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))-Fx*np.sin(alpha+g.alpha_i+g.alpha_0+g.ip)/(g.m*V*np.cos(beta))
    A[3:6]=np.dot(inv(I), np.array([Mt[0],Mt[1],Mt[2]])+F[3:6]-np.cross(np.array([p,q,r]),np.dot(I,np.array([p,q,r]))))
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7]=q*math.cos(phi) -r*math.sin(phi)
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)




    return A









def Constraints_Beta(x, fix, CoefMatrix, atmo, g, PropWing):

    """function defining constraints for power minimization
       DeltaR is a fixed parameter
       beta is part of the state vector x
       inputs:
        -x =[alpha, beta, p, q, r, phi, theta, delta_a, delta_e, delta_i]
        x is the state to determine
        length of x except the propulsion levels is 8
        -fix = [V, deltaR, gamma, omega]
        fix is the vector of parameters whom are fixed by the user
    """

    rho = atmo[1]
    # First thing to do is to determine the number of engines on each semi wing
    n_eng=int(g.N_eng/2)

    # --- Now prepare variables for equations ---
    V=fix[0]
    alpha=x[0]
    beta = x[1]
    deltaR = fix[1]
    gamma = fix[2]                                     # Correct, something wrong here, two defined with x[1]...
    omega = fix[-1]
    p = x[1]
    q = x[2]
    r = x[3]
    phi = x[4]
    theta = x[5]
    I = np.array([ [g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz] ])
    dx = np.copy(x[-g.N_eng:])
    
    # --- Compute aerodynamic forces ---
    #here subvector  must be : (alpha, beta, p, q, r, da, dr, de, dx)
    sub_vect=np.array([alpha,beta,p,q,r])
    if g.nofin==False:
        sub_vect=np.append(sub_vect,[x[6],deltaR,x[7]]) # rudder is allowed
    else:
        sub_vect=np.append(sub_vect,[x[6],0,x[7]]) # no fin allowed
    

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), dx, atmo, g, PropWing)     #  Check because PattersonAugmented works with Tc
    print("Forces :")
    print(F)

    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng

    #Now compute thrust and moments
    Fx_vec=g.Thrust(x[-g.N_eng:],V_vect)
    Fx = np.sum(Fx_vec)


    Moment= np.zeros((g.N_eng,3))
    for i in range(g.N_eng):
        a= np.array([ g.x_m , g.PosiEng[i] , g.z_m])
        b=np.array([ Fx_vec[i]*np.cos(g.alpha_i + g.alpha_0+g.ip)  ,  0  ,  -Fx_vec[i]*np.sin(g.alpha_i + g.alpha_0+g.ip)  ])
        Moment[i,:] = np.cross(a,b)
    Thrust_moment_body_axis =np.array(( np.sum(Moment[:,0]), np.sum(Moment[:,1]) , np.sum(Moment[:,2]) ) )

    Body2Aero_matrix = np.array([   [np.cos(alpha)*np.cos(beta), np.sin(beta) , np.sin(alpha)*np.cos(beta) ], [ -np.cos(alpha)*np.sin(beta) , np.cos(beta) , -np.sin(beta)*np.sin(beta) ] , [ -np.sin(alpha), 0   , np.cos(alpha)  ]])

    Trust_moment_aero_axis = Body2Aero_matrix @  Thrust_moment_body_axis

    Mt = Trust_moment_aero_axis
    


#     Now sum up the constraints:
    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi) 
    
    A=np.zeros(10+g.inop)
    """
    A0 = x
    A1 = y
    A2 = z
    A3 = l
    A4 = m
    A5 = n
    A6 = phi
    A7 = theta
    A8 = gamma
    A9 = Omega
    """
    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx*np.cos(alpha+g.alpha_i+g.alpha_0+g.ip)*np.cos(beta)/g.m
    A[1]=(p*np.sin(alpha) - r * np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m * V)
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))-Fx*np.sin(alpha+g.alpha_i+g.alpha_0+g.ip)/(g.m*V*np.cos(beta))
    A[3:6]=np.dot(inv(I), np.array([Mt[0],Mt[1],Mt[2]])+F[3:6]-np.cross(np.array([p,q,r]),np.dot(I,np.array([p,q,r]))))
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7]=q*math.cos(phi) - r * math.sin(phi)
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)

    for i in range(g.inop):
        A[-1-i]=x[-1-i]

    if g.hangar['version']=='original':
        #no DEP with original twin or N engines; all engines have the same thrust
        D=np.copy(A)
        for i in range(g.N_eng-g.inop-1):
            AAd=x[-g.N_eng]-x[-g.N_eng+i+1]
            D=np.append(D,[AAd])
        return D
    else:
        return A



def fobjective(x, fix, rho, g):
    
#Power=np.sum(x[-g.N_eng:])*2*g.P_var/float(g.N_eng)*rho/1.225/1000000
    Power=np.sum(x[-g.N_eng:])*2*g.P_var/float(g.N_eng)/1000000
    
    return Power



def fobjectiveBeta(x, fix, CoefMatrix, rho, g):
    """inputs:
        x =[alpha, beta, p, q, r, phi, theta, delta_a, delta_e, delta_i]"""
    Dx = x[-g.N_eng:]
    MeanDx = np.mean(Dx)
    
    BetaSideCoef = np.dot(CoefMatrix[1],x[0:7]) #side force
    
    return abs(BetaSideCoef)



def fobjectivePropWingInterac(x, fix, rho, g):

    Dx = x[-g.N_eng:]
    MeanDx = np.mean(Dx)
    stdDx = np.std(Dx)

    if g.hangar['aircraft']=='ATR72':
         #    Power=np.sum(x[-g.N_eng:])*2*g.P_var/float(g.N_eng)*rho/1.225/1000000
         return MeanDx+stdDx

    elif g.hangar['aircraft']=='DECOL':
         return MeanDx*0.5+stdDx*0.5



    
def fobjectiveBaseThetas(x,fix,rho,g):


    if g.hangar['aircraft']=='ATR72':
         #    Dthetas=x[-int(g.N_eng/2):]
         #    Power = np.dot(g.BasisErf(g.BasisCenters),np.concatenate((np.flip(Dthetas,axis=0),Dthetas)))*g.P_a/(g.b/2)/g.NBasisCenter
         Power = g.BasisFunPowerTotal(x[-int(g.N_eng/2):])
         #    stdThetas = np.std(Dthetas)
         #    VTsize = x[-int(g.N_eng/2)-1]
         #    dx = x[-int(g.N_eng/2)*3-1:-int(g.N_eng/2)-1]
         #    DeflRudder = 0#x[8]
    
         #    DeltaPVT = (VTsize-1)*0.107*1e6/(2*g.P_var)
    
         #    return np.max(Dthetas)
         #    return (np.mean(Dthetas)+stdThetas)/2
         return (Power/(g.P_a*2))

    elif g.hangar['aircraft']=='DECOL':

         Dthetas=x[-int(g.N_eng/2):]
         #    stdThetas = np.std(Dthetas)
         VTsize = x[-int(g.N_eng/2)-1]
         dx = [0,0]#x[-int(g.N_eng/2)*3-1:-int(g.N_eng/2)-1]
         DeflRudder = 0#x[8]

         return (sum(Dthetas)/(g.N_eng/2)+VTsize+sum(dx)/g.N_eng+DeflRudder)/3




def fobjective_minimum_alpha(x, fix, rho, g):

    alpha = x[0]

    return alpha



def fobjective_minimum_gamma(x, fix, rho, g):

    abs_gamma = abs(x[-2])


    return abs_gamma











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






















