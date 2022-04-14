
"""
Useful model for giving out the aerodynamic forces whichever the variables from state vector and fix vector are, so not
need for trim. If Aeroforces is changed it should be changed here too.

author: david.planas-andres

"""

import numpy as np
import math
from numpy.linalg import inv

import ReadFileUtils as Read  # utils to read Xfoil file
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import matplotlib.pyplot as plt


def Constraints_DEP(CoefMatrix, atmo, g, PropWing):



    #x = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i(hay 12), V , beta , gamma, omega] = x + fix


    rho = atmo[1]
    n_eng = int(g.N_eng / 2)
    PW = PropWing


    # --- Now prepare variables for equations ---
    V = 72
    alpha = 0.07891
    beta = 0
    gamma = 0
    omega = 0
    p = 0
    q = 0
    r = 0
    phi = -0.0000493944669
    theta = alpha+gamma
    aileron = 0.00006
    elevator = 0.01572
    rudder = 0.00022
    delta_x = 0.33001

    # Flaps : Flaps must be changed from main.

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




    Tc = Fx_vec / (2 * rho * g.Sp * V_vect ** 2)



    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi)

    fixtest = np.array([V, beta, gamma, omega])

    g.IsPropWing = True
    g.IsPropWingDrag = True


    F = CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    # F contiene solo las fuerzas y momentos aerodinámicos en ejes viento.


    printx(x, fixtest, atmo, g, PW)



    g.IsPropWing = False
    g.IsPropWingDrag = False

    F_no_int = CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)



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




    '''
    Coefs=np.zeros(len(F))

    for i in range(len(F)):
        if i==0 or i==1 or i==2:
            Coefs[i] = F[i] / (0.5 * rho * V**2 * g.S)

        elif i==4:
            Coefs[i] = F[i] / (0.5 * rho * V**2 * g.S * g.c)

        else:
            Coefs[i] = F[i] / (0.5 * rho * V ** 2 * g.S *  g.b)
    '''


    return F,F_no_int




def CalcForce_aeroframe_DEP(V, CoefMatrix, x, Tc, atmo, g, PropWing):
    """ Function to compute aerodynamic forces in the velocity frame (aero frame)
    for the DEP configuration. The propulsion force and moments are not computed here
    Since V is fixed, Coef Matrix must be calculated before
    Can handle DEP, with and without rudder, 2 or more engines
    """
    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]
    r = x[4]

    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)
    # set non dim for p,q,r
    nonDim = np.ones(len(x))
    nonDim[2] = g.b/(2*V)
    nonDim[3] = g.c/(2*V)
    nonDim[4] = g.b/(2*V)
    x = x*nonDim


    #If prop wing, aileron contributions are included in patterson; Cl_delta_a ; Cm_delta_a ; Cn_delta_a
    if g.IsPropWing and g.isail:
        CoefMatrix[3:6,5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0

    #    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F = np.zeros(3)
    M = np.zeros(3)
    xsym = np.copy(x)


    xsym[1] = abs(xsym[1])  # make beta always positive since derivatives have already correct sign for drag and lift only
    xsym[5] = abs(xsym[5])  # make ailerons deflection always positive for drag increase and lift decrease

    if g.nofin == False:
        xsym[-1] = abs(xsym[-1]) # make rudder deflection always positive for drag increase and lift decrease


    F[0] = np.dot(CoefMatrix[0, 1:], xsym[1:])               # (                CD_BETA*|BETA| + CD_P*P^ +  CD_Q*Q^ + CD_R*R^ + CD_DA*|DA| +  CD_DE*DE   + CD_DR*|DR|)  !not alpha
    F[1] = np.dot(CoefMatrix[1], x)                          # ( CY_ALFA*ALFA + CY_BETA* BETA  + CY_P*P^ +  CY_Q*Q^ + CY_R*R^ + CY_DA*DA   +  CY_DE*DE   + CY_DR*DR)
    F[2] = np.dot(CoefMatrix[2], xsym)                       # ( CL_ALFA*ALFA + CL_BETA*|BETA| + CL_P*P^ +  CL_Q*Q^ + CL_R*R^ + CL_DA*|DA| +  CL_DE*DE + CD_DR*|DR|)  !calculated again later if interaction
    M = np.dot(CoefMatrix[3:6, :], x)

    DragQuad = F[0] + g.Cda*x[0]**2 + g.Cdb * x[0] + g.Cdc

    g.Cm_alpha_aero_interaction = Cm_alpha(V, CoefMatrix, x, Tc, atmo, g, PropWing)

    if g.IsPropWing:

        CoefMatrix[4,0] = g.Cm_alpha_aero_interaction                                                                   #Use new CM_ALPHA
        M = np.dot(CoefMatrix[3:6,:],x)                                                                                 #Redo calculus with new Cm_alpha
        F[2] = np.dot(CoefMatrix[2][1:],xsym[1:]) + g.aht*x[0]-g.CL0_HT                                                 #F[2] Calculated without alpha: CL_BETA*BETA + CL_P*P + CL_Q*Q + CL_R*R + CL_DA*|DA| + CL_DE*DE + CL_DR*|DR|   )
        #Terms for horizontal tail added (alpha, and 0-alpha term) to modify x[0] and g.CL0_HT to take into account slisptream
        if V <= g.VelFlap:

            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g, beta, p, V, r)


            if len(CLCl)>2 and g.IsPropWingDrag:
                # Drag is computed by patterson, add contribution of other variables (than alpha and dx)
                Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]]) # add alpha=0 coefficients  #Need to revise if CD0T must be here or is counted in Patterson
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])


            else:
                Fbody = np.array([-DragQuad,F[1],-F[2]-CLCl[0]])  # add alpha=0 coefficients
                # add roll effect
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])

        else:
            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, 0, g, beta, p, V, r)
            # by default lift and drag are computed here
            Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]]) # add alpha=0 coefficients      #Need to revise if CD0T must be here or is counted in Patterson
            # add roll effect
            Moment = M+np.array([CLCl[1], g.Cm0, CLCl[4]])
    else:

        if V <= g.VelFlap:
            Fbody = np.array([-DragQuad, F[1], -F[2]-g.CL0 - g.CL0_fl])
            Moment = M+np.array([0, g.Cm0 + g.Cm0_fl, 0])
        else:
            Fbody = np.array([-DragQuad, F[1], -F[2]-g.CL0])  # add alpha=0 coefficients
            Moment = M+np.array([0, g.Cm0, 0])

    g.TotalDrag = abs(Fbody[0])
    g.Lift = abs(Fbody[2])
    g.lift2drag = abs(Fbody[2]/Fbody[0])
    Fbody = 0.5*V**2.0*rho*g.S*Fbody
    Moment = 0.5*V**2.0*rho*g.S*Moment*np.array([g.b,g.c,g.b])

    return np.append(Fbody, Moment)
















def Cm_alpha(V, CoefMatrix, x, Tc, atmo, g, PropWing):
    """ Function to compute aerodynamic forces in the velocity frame (aero frame)
    for the DEP configuration. The propulsion force and moments are not computed here
    Since V is fixed, Coef Matrix must be calculated before
    Can handle DEP, with and without rudder, 2 or more engines"""

    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]
    r = x[4]

    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)
    # set non dim for p,q,r
    nonDim = np.ones(len(x))
    nonDim[2] = g.b/(2*V)
    nonDim[3] = g.c/(2*V)
    nonDim[4] = g.b/(2*V)
    x = x*nonDim


    if g.IsPropWing and g.isail:
        CoefMatrix[3:6,5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0


    #Additional Parameters.



    g.Vef2Vinf_2 = PropWing.Augmented_velocity_wing(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g,beta,p,V,r)          #(V_ef/V_inf)^2  (is not tail/Vinf , keep that in mind)

    g.eta_downwash = 0.8
    # g.X_CA_wb = g.x_cg/g.c - (      (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *0.8)/ (CoefMatrix[2,0]-g.aht)   )   #--> Calculated so that downwash *ratio of dynamic pressures is 0.8



    g.SM_CA_wingfuselage =  (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *g.eta_downwash)/ (CoefMatrix[2,0]-g.aht)      #it bothers me that if you calculate the CA_wingbody does not really match with geometry,
    #maybe because center of gravity is not well stimated? IN ATR for having a good value for
    #CA_wingbody, center of gravity should be around 11.4 m, is 12.41 now. Anyway calculus of
    # Cm alpha only cares about static margin


    #Cl_alpha with interaction calculus
    alpha_1 = 0 * np.pi/180
    alpha_2 = 2 * np.pi/180
    CL_alpha_interaction = (   (PropWing.CalcCoef(Tc,V/a_sound, atmo, alpha_2,dail,g.FlapDefl,g,beta,p,V,r)[0] + g.aht*alpha_2 ) - (PropWing.CalcCoef(Tc,V/a_sound, atmo, alpha_1,dail,g.FlapDefl,g,beta,p,V,r)[0] + g.aht*alpha_1 ) )/ (alpha_2-alpha_1)


    #CALCULUS OF NEW CM_ALPHA_AERODYNAMIC

    Cm_alpha_interaction = (1 + ((CL_alpha_interaction - CoefMatrix[2,0]) * g.SM_CA_wingfuselage)/CoefMatrix[4,0]) * CoefMatrix[4,0]
    #this formula is valis supossing that downwash, c.gravity and tail-wing pressure ratio does not change when implementing DEP


    return Cm_alpha_interaction






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






















