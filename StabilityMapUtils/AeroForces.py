# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 10:51:13 2017

Module for coef interpolation and forces calculation

@author: e.nguyen-van
         david.planas-andres
"""
import numpy as np


def Interpol(V, A1, A2, v1, v2):
    # Function to interpol any kind of variables in function of the velocity V
    # input : 
    # V : current velocity
    # A1, A2 : lower matrix, higher matrix
    # v1, v2 : velocities corresponding to matrices
    a = (A2-A1)/(v2-v1)
    b = A1-a*v1
    Areturn = a*V+b
    return Areturn

def InterpolRho(V, rho, v):
    # function to interpol Rho, since altitude is function of Velocity
    if V < v[0]:
        return rho[0]
    
    elif V > v[-1]:
        return rho[-1]
    else:
        exitcondition = 1
        length_v = len(v)-1
        i = 0
        while exitcondition:
           
            if V == v[i]:
                rhoreturn = rho[i]
                exitcondition = 0
            
            elif V > v[i] and V < v[i+1]:
                rhoreturn = Interpol(V, rho[i], rho[i+1], v[i], v[i+1])
                exitcondition = 0  # exit
            
            else:
                i = i+1
                
            if i == length_v:  # security to exit the while
                print("AeroForces : Error in interpolating rho, returning 0")
                rhoreturn = 0
                exitcondition = 0
    
    return rhoreturn

def CoefInterpol( V, A, v):
    # A, function takes a numpy array composed of all matrices A [A1; A2; ...], all array types!!
    # v, an array of corresponding velocity
    # V, the velocity at which to compute the coef
    # size of each matrix 6 row, m column
    row = 6
    nmatrix = len(A[:, 1])/6
    if nmatrix == float:
        # error
        print('Aero forces: general coef matrix not a multiple of 6')
        return 0
    
    elif V < v[0]:
        # ill posed problem, the velocity is below smallest velocity for coef
        # use first matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is below first ref velocity for coef v = {1:0.2f}".format(V, v[0]))
        print("Continue with first matrix")
        return A[0:row, :]
    
    elif V > v[-1]:
        # same idea, the velocity is greater than max used to determine coef. Results in flight faster than cruise
        # use last matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is higher than last ref velocity for coef v = {1:0.2f}".format(V, v[-1]))
        print("Continue with last matrix")
        return A[-6:]


        
    else : #otherwise interpolate
        exitcondition = 1
        length_v = len(v)-1
        i = 0
        while exitcondition:
           
            if V == v[i]:
                Areturn = A[i*row:i*row+row]
                exitcondition = 0
            
            elif V > v[i] and V < v[i+1]:
                Areturn = Interpol(V, A[i*row:i*row+row, :], A[(i+1)*row:(i+1)*row+row, :], v[i], v[i+1])
                exitcondition = 0  # exit
            
            else:
                i = i+1
                
            if i == length_v+1:  # security to exit the while
                print("!!! FAILURE !!! AeroForces : Error in interpolation, returning 0")
                Areturn = 0
                exitcondition = 0

    return Areturn










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

    # Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, de, dr) 
    # set non dim for p,q,r
    nonDim = np.ones(len(x))
    nonDim[2] = g.b/(2*V)
    nonDim[3] = g.c/(2*V)
    nonDim[4] = g.b/(2*V)
    x = x*nonDim

    # If prop wing, aileron contributions are included in patterson; Cl_delta_a ; Cm_delta_a ; Cn_delta_a
    if g.IsPropWing and g.isail:
        CoefMatrix[3:6, 5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0

    F = np.zeros(3)
    M = np.zeros(3)

    xsym = np.copy(x)
    xsym[1] = abs(xsym[1])  # make beta always positive since derivatives have already correct sign for drag and lift only
    xsym[5] = abs(xsym[5])  # make ailerons deflection always positive for drag increase and lift decrease

    if g.nofin == False:
              xsym[-1] = abs(xsym[-1])  # make rudder deflection always positive for drag increase and lift decrease

    F[0] = np.dot(CoefMatrix[0, 1:], xsym[1:])               # (                CD_BETA*|BETA| + CD_P*P^ +  CD_Q*Q^ + CD_R*R^ + CD_DA*|DA| +  CD_DE*DE   + CD_DR*|DR|)  !not alpha
    F[1] = np.dot(CoefMatrix[1], x)                          # ( CY_ALFA*ALFA + CY_BETA* BETA  + CY_P*P^ +  CY_Q*Q^ + CY_R*R^ + CY_DA*DA   +  CY_DE*DE   + CY_DR*DR)
    F[2] = np.dot(CoefMatrix[2], xsym)                       # ( CL_ALFA*ALFA + CL_BETA*|BETA| + CL_P*P^ +  CL_Q*Q^ + CL_R*R^ + CL_DA*|DA| +  CL_DE*DE + CD_DR*|DR|)  !calculated again later if interaction
    M = np.dot(CoefMatrix[3:6, :], x)

    DragQuad = F[0] + g.Cda*x[0]**2 + g.Cdb * x[0] + g.Cdc + (g.CD0T - g.Cdc_fl_0)    #  Last term for moving above the polar Cd0. (CD0T more accurate than Cdc_fl_0)

    if g.IsPropWing:
        # Taking from Cl_beta, Cl_p, Cl_r the contribution of the wing, since effects from beta,p,r are re-calculated on wing's LIFT inside Patterson. Lift creates roll!
        CoefMatrix[3, 1] = CoefMatrix[3, 1] - g.Matrix_no_tail_terms[3, 1]
        CoefMatrix[3, 2] = CoefMatrix[3, 2] - g.Matrix_no_tail_terms[3, 2]
        CoefMatrix[3, 4] = CoefMatrix[3, 4] - g.Matrix_no_tail_terms[3, 4]
        M = np.dot(CoefMatrix[3:6, :], x)

        Cm, CL_tail, eps = Cm_and_CL_tail(V, CoefMatrix, x, Tc, atmo, g, PropWing)       #  For the pitching moment calculus with Delft paper
        F[2] = np.dot(CoefMatrix[2][1:], xsym[1:]) + CL_tail                                                # F[2] Calculated without alpha: CL_BETA*BETA + CL_P*P + CL_Q*Q + CL_R*R + CL_DA*|DA| + CL_DE*DE + CL_DR*|DR|   )
                                                                                                            # Terms for horizontal tail added (alpha, and 0-alpha term) to modify x[0] and g.CL0_HT to take into account slisptream
        if V <= g.VelFlap or g.FlapDefl != 0:

            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g, beta, p, V, r)

            Cd_Jam = Jamessondrag(V, CoefMatrix, x, Tc, atmo, g, PropWing, CLCl[0]+CL_tail)

            if len(CLCl) > 2 and g.IsPropWingDrag:
                # Drag is computed by patterson, add contribution of other variables (than alpha and dx)
                Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]])  # add alpha=0 coefficients
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])

                if g.hangar['aircraft'] == 'X-57':
                    Fbody[0] = -Cd_Jam - F[0] -CLCl[3] - CLCl[5]

            else:
                Fbody = np.array([-DragQuad, F[1], -F[2]-CLCl[0]])  # add alpha=0 coefficients
                # add roll effect
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])



        else:
            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, 0, g, beta, p, V, r)
            # by default lift and drag are computed here
            Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]])  # add alpha=0 coefficients
            # add roll effect
            Moment = M+np.array([CLCl[1], g.Cm0, CLCl[4]])
            Cd_Jam = Jamessondrag(V, CoefMatrix, x, Tc, atmo, g, PropWing, CLCl[0]+CL_tail)
            if g.hangar['aircraft'] == 'X-57':
                Fbody[0] = -Cd_Jam - F[0] - CLCl[3] - CLCl[5]
    else:

        if V <= g.VelFlap or g.FlapDefl != 0:
            Fbody = np.array([-DragQuad, F[1], -F[2]-g.CL0 - g.CL0_fl])
            Moment = M+np.array([0, g.Cm0 + g.Cm0_fl, 0])
        else:
            Fbody = np.array([-DragQuad, F[1], -F[2]-g.CL0])  # add alpha=0 coefficients
            Moment = M+np.array([0, g.Cm0, 0])

    if g.IsPropWing:
        Moment[1] = Cm
    g.TotalDrag = abs(Fbody[0])
    g.Lift = abs(Fbody[2])
    g.lift2drag = abs(Fbody[2]/Fbody[0])
    Fbody = 0.5*V**2.0*rho*g.S*Fbody
    Moment = 0.5*V**2.0*rho*g.S*Moment*np.array([g.b, g.c, g.b])

    return np.append(Fbody, Moment)















def Cm_and_CL_tail(V, CoefMatrix, x, Tc, atmo, g, PropWing):

    """
    Function for computing the total pitch moment and the horizontal tail's lift while taking into account the
    effects of the slipstream. Based on Obert's theory, a correlation of wind tunnel data.
    More can be found in:
             Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof; DOI 10.2514/6.2017-0236

    Here x must be of the form (alpha, beta, p, q, r, da, de, dr)

    We need to compute several special magnitudes:

          * VarVtoV : Relation between the augmentation of speed in the slipstream and free stream velocity. (The speed
                      in the slipstream is (V + VarV), so we do (V + VarV)/V )

          * D_s : slipstream tube diameter.

          * CL0w: (Wing + fus) lift at 0 angle of attack and 0 deflection of flaps in presence of slipstream.

          * alpha_w_0 : This is the alpha for which the lift of the wing is 0 when there are no flaps and no thrust (no slipstream).
                     Keep in mind this angle is not modified by having slipstream. Slipstream only modifies the slope of the curve,
                     but the all the different CL curves for different Ct start at this point (alpha_w_0,0). It is modified
                     if we deployed the flaps.

          * VarCLs0 :  This is the difference between
                  1) The wing lift at alpha_w_0  for the given flaps with slipstream (Ct =! 0)
                  2) The wing lift at alpha_w_0 withtout slipstream (Ct = 0). This is 0 if flaps retracted,
                     if flaps deployed this is exactly CL0_flaps

          * VarCLsalpha : This is the difference between:
                  1) The wing lift at the given alpha  for the given flaps with slipstream (Ct =! 0)
                  2) The wing lift at the given alpha for the given flaps withtout slipstream (Ct = 0).

    To extend, see Fig:20 in "The effect of propeller slipstream on the static longitudinal stability and control
    of multi-engined propeller aircraft" Ed. Obert


    @david_planas
    """

    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]/(g.b/(2*V))
    q = x[3]/(g.c/(2*V))
    r = x[4]/(g.b/(2*V))
    alpha = x[0]
    CL_alpha_no_int = CoefMatrix[2, 0]
    da = x[5]
    de = x[6]
    dr = x[7]

    Fx_vec = Tc * (2*rho*g.Sp*V**2)
    Fx = np.sum(Fx_vec)

    if g.IsPropWing and g.isail:
        CoefMatrix[3:6, 5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0


    # Establish how many engines in the area of influence of horizontal tail (are within the horizontal tail wingspan)
    engines = []
    for i in range(len(g.yp)):
        if abs(g.yp[i]) < (g.bh/2):
            engines.append(i)
    engines = np.array(engines)


    # Slipstream velocity to free stream velocity of each engine (vector)
    VarVtoV = (1+Fx_vec/(0.5*rho*g.Sp*V**2))**0.5 - 1

    # Contracted slipstream diameter of each engine (vector)
    D_s = g.Dp * ((V + 0.5 * V * VarVtoV)/(V + V * VarVtoV)) ** 0.5


    # (Wing + fus) lift at alpha = 0 and DeflFlaps = 0 in presence of slipstream.
    CL0w = PropWing.CalcCoef(Tc, V/a_sound, atmo, 0, dail, 0, g, beta, p, V, r)[0]

    # Calculus of Wing zero-lift angle of attack (Is the same with or without slipstream), no flaps
    alpha_0_w = -(g.CL0-g.CL0_HT)/(CL_alpha_no_int - g.aht)  # [rad]

    # At alpha_w_0, for the wing, for the given FlapDefl, the difference of lift between the case with interaction and without it (Ct =!0 and Ct =0)
    VarCLs0 = (PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha_0_w, dail, g.FlapDefl, g, beta, p, V, r)[0]) - ((CL_alpha_no_int - g.aht) * alpha_0_w + (g.CL0-g.CL0_HT) + g.CL0_fl)

    # For the wing, for the given FlapDefl, for the given alpha, the difference of lift between the case with interaction and without it (Ct =!0 and Ct =0)
    VarCLsalpha = (PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha, dail, g.FlapDefl, g, beta, p, V, r)[0])-((CL_alpha_no_int - g.aht) * alpha + (g.CL0-g.CL0_HT) + g.CL0_fl)

    # CL wing at the given alpha and FlapDefl condition, with slipstream
    CL = (PropWing.CalcCoef(Tc, V/a_sound, atmo, alpha, dail, g.FlapDefl, g, beta, p, V, r)[0])

    # Computing Tail's lift coefficient and Tail's pitching moment
    CL_tail, Cm_tail, eps = VerticalTail_Lift_and_moment(V, alpha, g, PropWing, CL, CL_alpha_no_int, D_s, VarVtoV, engines)

    # Computing tail-off pitching moment
    Cm_tail_off = Tail_off_Pitching_Moment(x, CoefMatrix, V, alpha, g, PropWing, CL_alpha_no_int, CL0w, D_s, VarVtoV, CL_tail, VarCLs0, VarCLsalpha)

    Cm = Cm_tail_off + Cm_tail

    return Cm, CL_tail, eps











def VerticalTail_Lift_and_moment(V, alpha, g, PropWing, CL, CL_alpha_no_int, D_s, VarVtoV, engines):

    """
    Function to compute the tail lift and pitching moment. The function computes the slipstream and the dynamic pressure
    in the horizontal tail when there is a propeller in front of the wing.

    The normal downwash would be:
    eps = g.eps0 + g.deps_dalpha * alpha

    With slipstream it would be:

    eps = (g.deps_dalpha / (CL_alpha_no_int-g.aht)) * CL
    """



    # We are just interested on contracted slipstream diameter and slipstream speed of the engines in the area of influence of HT
    VarVtoV = np.mean(VarVtoV[engines[0]:engines[-1]])  # Not anymore a vector, but a float
    D_s = np.mean(D_s[engines[0]:engines[-1]])  # Not anymore a vector, but a float


    if g.hangar['aircraft'] == 'X-57':
       eps_noinflow = (g.deps_dalpha / (CL_alpha_no_int-g.aht)) * (CL) + 0.88*np.pi/180 # For taking into account downwash at zero lift. Maybe equation eps0 + deps/dalpha*alpha was more precise

    else:
       # Epsilon without inflow effects
       eps_noinflow = (g.deps_dalpha / (CL_alpha_no_int-g.aht)) * (CL)

    # K_epsilon calculus
    if ((g.lh2/g.c) < 5):
        g.K_e = -0.0085*(g.lh2/g.c)**3 + 0.1078*(g.lh2/g.c)**2 - 0.5579*(g.lh2/g.c) + 2.4546
    else:
        g.K_e = 1.3



    # H calculus :   Vertical Distance of Slipstream Center Line to Horizontal Tail
    h = g.z_h_w - g.lh*np.sin(alpha) - (g.x_offset[int(g.N_eng/2)] + 0.25 * g.c)*np.sin(alpha) + g.lh2*np.sin(g.K_e * eps_noinflow) + g.FlChord * g.c * np.sin(g.FlapDefl) + 0.25 * (g.x_offset[int(g.N_eng/2)] + 0.25 * g.c) * np.sin(PropWing.alpha0_fl * g.FlapDefl)

    # PropWing.alpha0_fl is the change of alpha_0 for unitary deflection of flap. In °/° = rad/rad
    # If you want radians, it is necessary to multiply by flaps deflection, in radians.
    # g.FlapDefl is in radians

    if (h/(0.5*D_s))< 2.5 and (h/(0.5*D_s)) > -1:
        var_eps = -0.2263 * (h/(0.5*D_s)) ** 6 + 1.0584 * (h/(0.5*D_s)) ** 5 - 0.2971 * (h/(0.5*D_s)) ** 4 - 3.56 * (h/(0.5*D_s)) ** 3 + 0.7938 * (h/(0.5*D_s)) ** 2 + 5.6374 * (h/(0.5*D_s)) + 0.0246
        extra_eps = (var_eps * VarVtoV) * np.pi / 180
        eps = eps_noinflow + extra_eps
    elif (h/(0.5*D_s)) > -2.5 and (h/(0.5*D_s)) < -1:
        var_eps = y = -2.2973*(h/(0.5*D_s))**6 - 25.197*(h/(0.5*D_s))**5 - 111.87*(h/(0.5*D_s))**4 - 255.44*(h/(0.5*D_s))**3 - 314.31*(h/(0.5*D_s))**2 - 198.75*(h/(0.5*D_s)) - 53.729
        extra_eps = (var_eps * VarVtoV) * np.pi / 180
        eps = eps_noinflow + extra_eps
    else:
        var_eps = 0
        eps = eps_noinflow


    # Dynamic pressure ratio in the horizontal tail
    if (1 - (2*h / D_s)**2) > 0:
        bs = D_s * (1 - (2*h / D_s)**2) ** 0.5
        Sh_s = len(engines) * bs * g.c_ht
        dpratio = ((Sh_s / g.Sh) * (1 + VarVtoV)**2 + (1-Sh_s/g.Sh))

    else:
        dpratio = 1



    if g.hangar['aircraft'] == 'X-57':

        if (alpha + g.it - eps) < (16*np.pi/180):
            # TAIL LIFT
            CL_tail = g.aht2 * (alpha + g.it - eps) * dpratio
            # TAIL MOMENT
            Cm_tail = -CL_tail * (g.lv)/(g.c)


        elif (alpha + g.it - eps)>= (16*np.pi/180) and (alpha + g.it - eps)<(18*np.pi/180):
            # TAIL LIFT
            CL_tail = g.aht2 * (15*np.pi/180) * dpratio
            # TAIL MOMENT
            Cm_tail = -CL_tail * (g.lv)/(g.c)

        else:
            # TAIL LIFT
            CL_tail = g.aht2 * dpratio * (2*(16*np.pi/180) - (alpha + g.it - eps))
            # TAIL MOMENT
            Cm_tail = -CL_tail * (g.lv)/(g.c)


    else:
        # TAIL LIFT
        CL_tail = g.aht2 * (alpha + g.it - eps) * dpratio
        # TAIL MOMENT
        Cm_tail = -CL_tail * (g.lv)/(g.c)





    # Basically, applying more Ct creates higher downwash, this means that for high angles of attack, the detachment
    # occurs later in the horizontal tail.


    return CL_tail, Cm_tail, eps










def Tail_off_Pitching_Moment(x, CoefMatrix, V, alpha, g, PropWing, CL_alpha_no_int, CL0w, D_s, VarVtoV, CL_tail, VarCLs0, VarCLsalpha):

    """
    Function to compute the tail-off pitching moment.
    """

    # Preparations

    if g.hangar['aircraft'] == 'X-57':
        c_flaps = 0.7794  # Flap Fowler, calculated manually for the 30° deflection
    else:
        # Computes augmented chord when flaps are deflected, by Pithagoras
        c_flaps = g.c * np.sqrt(((1-g.FlChord) + g.FlChord*np.cos(g.FlapDefl))**2 + (g.FlChord*np.sin(g.FlapDefl))**2)




    # Tail-off clean pitching moment
    Cm1 = g.Cm0_wo_HT + g.Cm0_fl + g.Cm_alpha_wb*alpha + np.dot(CoefMatrix[4, 1:8], x[1:8])


    # 1st contribution. Moment due to augmented dynamic pressure on the profile
    Cm2 = ((D_s * g.c)/g.S) * g.cm_0_s * ((1+VarVtoV) ** 2 - 1)
    Cm2 = sum(Cm2)

    # 2nd contribution. Change in pitching moment when deploying flaps
    Cm3 = (c_flaps/g.c)*(-0.25 + 0.32 * (g.FlChord*g.c / c_flaps)) * (1+0.2*(1-np.sqrt(2) * np.sin(g.FlapDefl))) * VarCLs0


    Var_xac_fus = -0.25

    # PAPER AIAA

    if g.FlapDefl == 0:
        F = 0
    elif g.FlapDefl <= 30*np.pi/180:
        F = (0.5 * (g.FlapDefl*180/np.pi) / 30 + 0.25)*(c_flaps/g.c-1) + 0.05 * (g.FlapDefl*180/np.pi) / 30 + (Var_xac_fus / g.c) * (1-(g.FlapDefl*180/np.pi) / 30)
    else:
        F = -0.75 * (c_flaps/g.c-1) - 0.05

    Cm4 = F * VarCLsalpha


    # THESIS AND OBERT

    Cm5 = -0.25*(c_flaps/g.c - 1)*VarCLsalpha

    Cm6 = -(0.05 + 0.5*(c_flaps/g.c - 1))*VarCLsalpha


    # To take the moments in the centre of gravity, not in the aerodynamic point
    # Passing pitching moment due to wing + fuselage lift from aerdynamic center (g.lemac + 0.25*g.c) tp center of gravity (g.x_cg)
    # We just need to account for the extra wing lift due to slipstream!
    Cm7 = - (VarCLsalpha)*(g.lemac + 0.25*g.c - g.x_cg)/g.c


    # AIAA PAPER
    Cm_tail_off1 = Cm1 +Cm2 + Cm3 + Cm4 + Cm7

    # OBERT
    Cm_tail_off2 = Cm1 +Cm2 + Cm3 + Cm5 + Cm6 + Cm7

    # Differences are very small so we take Oberts method and we do not need to provide Var_xac_fus = -0.25


    return Cm_tail_off2







def Jamessondrag(V, CoefMatrix, x, Tc, atmo, g, PropWing, CL):
    """
    For computing drag in all kind of cases, wether there are flaps and/or blowing into the wing.
    Inputs: (V, CoefMatrix, x, Tc, atmo, g, PropWing, CL)
    """

    rho = atmo[1]
    Fx_vec = Tc * (2*rho*g.Sp*V**2)
    Fx = np.sum(Fx_vec)

    # Slipstream velocity to free stream velocity of each engine (vector)
    VarVtoV = (1+Fx_vec/(0.5*rho*g.Sp*V**2))**0.5 - 1

    # Contracted slipstream diameter of each engine (vector)
    D_s = g.Dp * ((V + 0.5 * V * VarVtoV)/(V + V * VarVtoV)) ** 0.5

    VarVtoV_av = (sum(D_s*VarVtoV) / (sum(D_s)))
    mu = V/(V + VarVtoV_av*V)

    ARmu = (g.b**2/g.S) * (1+g.N_eng*mu**2)/(g.N_eng + mu**2)

    k = (mu**2) / (np.pi*0.8*ARmu)

    CD0 = 0.0537249 + (0.0782 - 0.0537249)*(g.FlapDefl)/(30*np.pi/180)


    CD = CD0 + k*CL**2

    return CD


