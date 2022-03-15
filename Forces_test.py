
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



    # --- Now prepare variables for equations ---
    V = 72
    alpha = 0.0665775709
    beta = 0
    gamma = 0.13020895
    omega = 0
    p = 0
    q = 0
    r = 0
    phi = -0.0000493944669
    theta = alpha+gamma
    aileron = 0.00001
    elevator = -0.02005
    rudder = 0.00020
    delta_x = 1

    x = np.array([alpha, p, q, r, phi, theta, aileron , elevator , rudder])

    for i in range(int(g.N_eng)):
       x= np.append(x,delta_x )

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



    g.IsPropWing = True
    g.IsPropWingDrag = True


    F = CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    # F contiene solo las fuerzas y momentos aerodinámicos en ejes viento.


    g.IsPropWing = False
    g.IsPropWingDrag = False

    F_no_int = CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)


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
    beta=x[1]
    p=x[2]
    r=x[4]

    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)
    #Determine correction for local air velocity at engines

    # set non dim for p,q,r
    nonDim=np.ones(len(x))
    nonDim[2]=g.b/(2*V)
    nonDim[3]=g.c/(2*V)
    nonDim[4]=g.b/(2*V)
    x=x*nonDim


    #If prop wing, aileron contributions are included in patterson; Cl_delta_a ; Cm_delta_a ; Cn_delta_a
    if g.IsPropWing and g.isail:
        CoefMatrix[3:6,5] = np.zeros(3)
        dail = x[5]
    else:
        dail=0

    #    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F=np.zeros((3))
    M=np.zeros((3))
    xsym=np.copy(x)

    #
    xsym[1]=abs(xsym[1]) # make beta always positive since derivatives have already correct sign for drag and lift only
    xsym[5]=abs(xsym[5]) # make ailerons deflection always positive for drag increase and lift decrease

    if g.nofin==False:
        xsym[-1]=abs(xsym[-1]) # make rudder deflection always positive for drag increase and lift decrease



    F[0]=np.dot(CoefMatrix[0,1:],xsym[1:])                                                                              # (                CD_BETA*|BETA| + CD_P*P^ +  CD_Q*Q^ + CD_R*R^ + CD_DA*|DA| +  CD_DE*DE   + CD_DR*|DR|)  !NOT ALPHA
    F[1]=np.dot(CoefMatrix[1],x) #side force                                                                            # ( CY_ALFA*ALFA + CY_BETA* BETA  + CY_P*P^ +  CY_Q*Q^ + CY_R*R^ + CY_DA*DA   +  CY_DE*DE   + CY_DR*DR)
    F[2]=np.dot(CoefMatrix[2],xsym)                                                                                     # ( CL_ALFA*ALFA + CL_BETA*|BETA| + CL_P*P^ +  CL_Q*Q^ + CL_R*R^ + CL_DA*|DA| +  CL_DE*|DE| + CD_DR*|DR|)  !CALCULATED LATER AGAIN IF INTERACTION
    M=np.dot(CoefMatrix[3:6,:],x)

    #Drag force quadratic in alpha
    #    DragQuad = F[0]+0.8029*x[0]**2 + x[0] * 0.12537 + 0.006434
    #     drag force from identified DECOL:
    #    DragQuad = F[0] + 0.99747*x[0]**2 + 0.18081*x[0] +0.02977
    DragQuad = F[0] + g.Cda*x[0]**2 + x[0] * g.Cdb + g.Cdc        # g.Cda = 1.5660758    g.Cdb = 0.03718442    g.Cdc = 0.00840314    ADD the Cd_alpha*alpha but with a polar.









    #Additional Parameters.


    g.Hor_tail_coef_vol = (g.Sh*g.lv) / (g.S*g.c)                                                                       #1.02   In a master thesis from Hamburg he gets 1.05 so I validate the value.
    g.eta_downwash = 0.8
    # g.X_CA_wb = g.x_cg/g.c - (      (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *0.8)/ (CoefMatrix[2,0]-g.aht)                         )    #--> Calculated so that downwash *ratio of dynamic pressures is 0.8

    #it bothers me that if you calculate the CA_wingbody does not really match with geometry,
    #maybe because center of gravity is not well stimated? IN ATR for having a good value for
    #CA_wingbody, center of gravity should be around 11.4 m, is 12.41 now. Anyway calculus of
    # Cm alpha only cares about static margin



    g.SM_CA_wingfuselage =  (CoefMatrix[4,0] + g.aht*g.Hor_tail_coef_vol *g.eta_downwash)/ (CoefMatrix[2,0]-g.aht)

    #Cl_alpha with interaction calculus
    alpha_1 = 0 * np.pi/180
    alpha_2 = 2 * np.pi/180
    CL_alpha_interaction = (   (PropWing.CalcCoef(Tc,V/a_sound, atmo, alpha_2,dail,g.FlapDefl,g,beta,p,V,r)[0] + g.aht*alpha_2 ) - (PropWing.CalcCoef(Tc,V/a_sound, atmo, alpha_1,dail,g.FlapDefl,g,beta,p,V,r)[0] + g.aht*alpha_1 ) )/ (alpha_2-alpha_1)


    #CALCULUS OF NEW CM_ALPHA_AERODYNAMIC

    g.Cm_alpha_aero_interaction =( 1 + ((CL_alpha_interaction - CoefMatrix[2,0]) * g.SM_CA_wingfuselage)/CoefMatrix[4,0] ) * CoefMatrix[4,0]
    #this formula is valis supossing that downwash, c.gravity and tail-wing pressure ratio
    # does not change when implementing DEP









    if g.IsPropWing:

        CoefMatrix[4,0]=g.Cm_alpha_aero_interaction
        M=np.dot(CoefMatrix[3:6,:],x)

        if V<=g.VelFlap:
            CLCl = PropWing.CalcCoef(Tc,V/a_sound, atmo, x[0],dail,g.FlapDefl,g,beta,p,V,r) #flap reduced to 15°
            # compute effects of other variables
            F[2]=np.dot(CoefMatrix[2][1:],xsym[1:])

            if len(CLCl)>2 and g.IsPropWingDrag:

                Fbody=np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T,F[1],-F[2]-CLCl[0]]) # add alpha=0 coefficients
                Moment=M+np.array([CLCl[1],g.Cm0+g.Cm0_fl*g.FlapDefl/15-CoefMatrix[4,0]*g.alpha_i,CLCl[4]])

            else:

                Fbody=np.array([-DragQuad-g.Cd0_fl-g.CD0T,F[1],-F[2]-CLCl[0]]) # add alpha=0 coefficients

                Moment=M+np.array([CLCl[1],g.Cm0_fl-CoefMatrix[4,0]*g.alpha_i,CLCl[4]])
                print("Yo")
        else :
            CLCl = PropWing.CalcCoef(Tc,V/a_sound, atmo, x[0],dail,0,g,beta,p,V,r)

            F[2]=np.dot(CoefMatrix[2][1:],xsym[1:])

            Fbody=np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T,F[1],-F[2]-CLCl[0]]) # add alpha=0 coefficients

            Moment=M+np.array([CLCl[1],g.Cm0-CoefMatrix[4,0]*g.alpha_i,CLCl[4]])
    else:
        if V<=g.VelFlap :
            DragQuad = F[0] + g.Cdafl*x[0]**2 + x[0] * g.Cdbfl + g.Cdcfl

            Fbody=np.array([-DragQuad-g.CD0T,F[1],-F[2]-g.CL0_fl-g.CL0])

            Moment=M+np.array([0,g.Cm0+g.Cm0_fl,0])
        else:
            Fbody=np.array([-DragQuad-g.CD0T,F[1],-F[2]-g.CL0])
            Moment=M+np.array([0,g.Cm0,0])


    #add contribution of horizontal tail

    if g.IsPropWing:
        Fbody[2]=Fbody[2]-g.aht*x[0]

    g.TotalDrag = abs(Fbody[0])
    g.Lift = abs(Fbody[2])
    g.lift2drag = abs(Fbody[2]/Fbody[0])
    Fbody=0.5*V**2.0*rho*g.S*Fbody
    Moment=0.5*V**2.0*rho*g.S*Moment*np.array([g.b,g.c,g.b])


    return np.append(Fbody, Moment)



