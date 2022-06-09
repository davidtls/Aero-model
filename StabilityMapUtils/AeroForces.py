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
        CoefMatrix[3:6, 5] = np.zeros(3)
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
              xsym[-1] = abs(xsym[-1])  # make rudder deflection always positive for drag increase and lift decrease


    F[0] = np.dot(CoefMatrix[0, 1:], xsym[1:])               # (                CD_BETA*|BETA| + CD_P*P^ +  CD_Q*Q^ + CD_R*R^ + CD_DA*|DA| +  CD_DE*DE   + CD_DR*|DR|)  !not alpha
    F[1] = np.dot(CoefMatrix[1], x)                          # ( CY_ALFA*ALFA + CY_BETA* BETA  + CY_P*P^ +  CY_Q*Q^ + CY_R*R^ + CY_DA*DA   +  CY_DE*DE   + CY_DR*DR)
    F[2] = np.dot(CoefMatrix[2], xsym)                       # ( CL_ALFA*ALFA + CL_BETA*|BETA| + CL_P*P^ +  CL_Q*Q^ + CL_R*R^ + CL_DA*|DA| +  CL_DE*DE + CD_DR*|DR|)  !calculated again later if interaction
    M = np.dot(CoefMatrix[3:6, :], x)

    DragQuad = F[0] + g.Cda*x[0]**2 + g.Cdb * x[0] + g.Cdc + (g.CD0T - g.Cdc_fl_0)    #  Last term for moving above the polar Cd0

    Cm, CL_tail = Cm_and_CL_tail(V, CoefMatrix, x, Tc, atmo, g, PropWing)       #  For the pitching moment calculus with Delft paper

    if g.IsPropWing:

        CoefMatrix[3, 2] = CoefMatrix[3, 2] - g.Matrix_no_tail_terms[3, 2]
        CoefMatrix[3, 4] = CoefMatrix[3, 4] - g.Matrix_no_tail_terms[3, 4]
        M = np.dot(CoefMatrix[3:6, :], x)

        F[2] = np.dot(CoefMatrix[2][1:], xsym[1:]) + CL_tail                                                #F[2] Calculated without alpha: CL_BETA*BETA + CL_P*P + CL_Q*Q + CL_R*R + CL_DA*|DA| + CL_DE*DE + CL_DR*|DR|   )
                                                                                                                        #Terms for horizontal tail added (alpha, and 0-alpha term) to modify x[0] and g.CL0_HT to take into account slisptream
        if V <= g.VelFlap or g.FlapDefl != 0:

            CLCl = PropWing.CalcCoef(Tc, V/a_sound, atmo, x[0], dail, g.FlapDefl, g, beta, p, V, r)


            if len(CLCl)>2 and g.IsPropWingDrag:
                # Drag is computed by patterson, add contribution of other variables (than alpha and dx)
                Fbody = np.array([-F[0]-CLCl[2]-CLCl[3]-CLCl[5]-g.CD0T, F[1], -F[2]-CLCl[0]]) # add alpha=0 coefficients
                Moment = M+np.array([CLCl[1], g.Cm0 + g.Cm0_fl, CLCl[4]])


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



    # From OpenVSP


    rho = atmo[1]
    a_sound = atmo[0]
    beta = x[1]
    p = x[2]
    q = x[3]
    r = x[4]
    alpha = x[0]
    CL_alpha_no_int = CoefMatrix[2, 0]
    da = x[5]
    de = x[6]
    dr = x[7]


    Fx_vec = Tc * (2*rho*g.Sp*V**2)
    Fx = np.sum(Fx_vec)


    # here x must be of the form (alpha, beta, p, q, r, da, de, dr)



    if g.IsPropWing and g.isail:
        CoefMatrix[3:6, 5] = np.zeros(3)
        dail = x[5]
    else:
        dail = 0




    #Additional Parameters.


    #Cl_alpha with interaction calculus
    CL_alpha_interaction = ((PropWing.CalcCoef(Tc, V/a_sound, atmo, 2 * np.pi/180, dail, g.FlapDefl, g, beta, p, V, r)[0] + g.aht*2*np.pi/180) - (PropWing.CalcCoef(Tc, V/a_sound, atmo, 0, dail, g.FlapDefl, g,beta, p, V, r)[0] + g.aht*0)) / ((2-0)*np.pi/180)
    CL0_int = (PropWing.CalcCoef(Tc, V/a_sound, atmo, 0, dail, 0, g,beta, p, V, r)[0] + g.aht*0 + g.CL0_HT)


    g.Vep2Vinf = (PropWing.Augmented_velocity_wing(Tc, V / a_sound, atmo, x[0], dail, g.FlapDefl, g, beta, p, V, r)) ** 0.5     # (V_ef/V_inf) (is not vtail/Vinf , keep that in mind)


    """
    g.eta_downwash = 0.858675
    g.X_CA_wb = g.x_cg/g.c - ((CoefMatrix[4, 0] + g.aht*g.Hor_tail_coef_vol * 0.8) / (CoefMatrix[2, 0]-g.aht))  
    g.SM_CA_wingfuselage = (CoefMatrix[4, 0] + g.aht*g.Hor_tail_coef_vol * g.eta_downwash) / (CoefMatrix[2, 0]-g.aht)

    #CALCULUS OF NEW CM_ALPHA_AERODYNAMIC

    Cm_alpha_interaction = (1 + ((CL_alpha_interaction - CoefMatrix[2, 0]) * g.SM_CA_wingfuselage)/CoefMatrix[4, 0]) * CoefMatrix[4, 0]
    #this formula is valid supossing that downwash, c.gravity and tail-wing pressure ratio does not change when implementing DEP
    """



    # Slipstream velocity to free stream velocity

    VarVtoV = (1+Fx/(0.5*g.N_eng*rho*g.Sp*V**2))**0.5 - 1  # Similar value than before with Patterson, is momenthum theory ...


    # Contracted slisptream diameter

    D_s = g.Dp * ((V + 0.5 * V * VarVtoV)/(V + V * VarVtoV)) ** 0.5


    # Epsilon without inflow effects

    eps = g.eps0 + (g.deps_dalpha / CL_alpha_no_int) * CL_alpha_interaction * alpha


    # H calculus :   Vertical Distance of Slipstream Center Line to Horizontal Tail

    h = g.z_h_w - g.lh*np.sin(alpha) - (g.xp + 0.25 * g.c)*np.sin(alpha) + g.lh2*np.sin(g.K_e * eps) + g.FlChord * g.c * np.sin(g.FlapDefl) + 0.25 * (g.xp + 0.25 * g.c) * np.sin(PropWing.alpha0_fl * g.FlapDefl)

    # PropWing.alpha0_fl has to be the change in airfoil section zero-lift angle-of-attack due to flap deflection (obtained using AVL)
    # PropWing.alpha0_fl is the change in alpha_0 for every radian deflection of flap. It is necessary to multiply by flaps deflection, in radians.
    # g.FlapDefl is in radians



    # Epsilon with inflow effects

    if h < 1.25 * D_s:
        extra_eps = (g.var_eps * VarVtoV) * np.pi / 180
        eps = eps + extra_eps



    # Dynamic pressure ratio in the horizontal tail

    V2 = (1 + VarVtoV) * V

    if (1 - (2*h / D_s)**2) > 0:
        bs = D_s * (1 - (2*h / D_s)**2) ** 0.5
        Sh_s = 2 * bs * g.c_ht

        dpratio = ((Sh_s / g.Sh) * (1 + VarVtoV)**2 + (1-Sh_s/g.Sh))
    else:
        dpratio = (1+VarVtoV)**2

    #BASICAMENTE PARECE QUE LA PARTE MOJADA DE LA COLA HORIZONTAL POR EL TUBO DE SLISPTREAM ES 0. ELLO IMPLICARIA QUE LA PRESION DINAMICA NO CAMBIA ?
    #MIRAR EN EL PAPER; VAYA MIERDA; EL SLIPSTREAM TAMPOCO CAMBIA O QUE...? EN DECOL SI AFECTARA QUAND MEME






    # TAIL-OFF PITCHING MOMENT!

    Cm_0 = g.Cm0_wo_HT + g.Cm0_fl + g.Cm_alpha_wb*alpha + np.dot(CoefMatrix[4, 1:8], x[1:8])




    c_flaps = g.c * np.sqrt(((1-g.FlChord) + g.FlChord*np.cos(g.FlapDefl))**2 + (g.FlChord*np.sin(g.FlapDefl))**2)

    Cm_s_0 = g.N_eng * ((D_s * g.c)/g.S) * g.cm_0_s * ((g.Vep2Vinf * V / V) ** 2 - 1)

    Cm_s_df = (c_flaps/g.c)*(-0.25+0.32*g.FlChord / c_flaps) * (1+0.2*(1-np.sqrt(2) * np.sin(g.FlapDefl))) * g.CL0_fl



    if g.FlapDefl == 0:

        F = 0

    elif g.FlapDefl <= 30*np.pi/180:

       F = (0.5 * (g.FlapDefl*180/np.pi) / 30 + 0.25)*(c_flaps/g.c-1) + 0.05 * (g.FlapDefl*180/np.pi) / 30 + (g.Var_xac_fus / g.c) * (1-(g.FlapDefl*180/np.pi) / 30)
    else:

       F = -0.75 * (c_flaps/g.c-1) - 0.05


    Cm_s_alpha = F * CL_alpha_interaction * alpha



    Cm_tail_off = Cm_s_0 + Cm_s_df + Cm_s_alpha + Cm_0 - (CL_alpha_interaction * alpha + g.CL0_fl + (CL0_int-g.CL0))*(g.lemac + 0.25*g.c - g.x_cg)



    # TAIL MOMENT
    Cm_tail = -(alpha + g.it - eps) * (g.aht2 * g.S/g.Sh) * dpratio * (g.Sh * g.lv)/(g.S * g.c)


    Cm = Cm_tail_off + Cm_tail



    # TAIL LIFT

    CL_tail = g.aht2 * (alpha + g.it - eps) * dpratio




    """ Function to compute longitudinal stability issues regarding the influence on the horizontal tail
    when there is interaction. To compute for the slipstream velocity and downwash and see the influence on
    the horizontal tail.

    Normally
    
    Cm_aero = Cm0 + Cm_alpha * alpha + Cm_q * q + Cm_de * de + Cm_u * u + Cm_d(alfa) * d(alfa) + Cm_d(de) * d(de)

    (Propulsion moments are already accounted in the equation modulus.)


                        -Cm_u = 0 as M is regime is always incompressible
                        -Cm_q is calculated in OpenVSP, more about how to calculate it in Page 345
                        -Cm_de          Page 345
                        -Cm_d(de)       Page 345
                        -Cm_d(alfa) ----> Teoría del retardo de estela ?




    Para proceder vamos a:

    Calcular el momento aerodinamico alrededor del centro aerodinamico conjunto ala-fuselaje. Este momento NO depende
    del angulo de ataque. Simulacion de OpenVSP con el ala y el fuselaje a 0 grados.


                    ATR: Without horizontal tail, 70m/s, Mach=0.2058 , alpha=0, ailerons,rudder, elevators ... = 0, XCG=11.571
    
                        Cm0 wingbody:  -0.171170
                        CL0 wingbody:  0.544738 this gives a moment
                        CD0 wingbody:  0.024993 this gives a moment


    Calcular el momento aerodinamico alrededor del centro aerodinamico de la cola horizontal. Este momento NO depende
    del angulo de ataque. Simulacion de OpenVSP con la cola horizontal a 0 grados.
    
    
                    ATR: Just horizontal tail
                    Cm0 = Cmo:                -0.001341
                    CL0 =  CLo:               -0.006814 This gives a moment, big one
                    CD0 =  CDo:                0.001353 a bit more this gives a moment


    Luego habra que calcular momentos de la sustentacion y la resistencia en el conjunto ala fuselaje
        Para ello:
                    Sustentacion en el cjto ala fuselaje ------>  (Lo tienes)

                    Resistencia en el cjto ala fuselaje -------> (lo tienes aunque hay que adaptarlo para que no
                    cuente solo la del aumento de presion dinamica en la parasita)

                    Centro aerodinamico en el cjto ala fuselaje. El fuselaje casi no sustenta, asi que se puede
                    considerar que sera el mismo que el del ala. Los tienes calculados, sino de todas formas OpenVSP deberia
                    dartelo

                    Calcular la sustentacion y la resistencia en el estabilizador horizontal.
                            Para ello : Necesitas un CL0 -----> de OpenVSP
                                        Necesitas un Cl_alfa -----> de OpenVSP
                                        Necesitas un angulo de ataque
                                        Necesitas una velocidad en la cola

                    Centro aerodinamico del estabilizador horizontal

                    Teoria para el calculo de la deflexion de estela

                    Teoria para el calculo de la presion dinamica en la cola.


    Y después calculamos directamente todo el momento aerodinamico EN EJES CUERPO Y LO CAMBIAMOS A EJES VELOCIDAD!! como:

    M = Cmac_wb + Cmac_HT + cg-ca ^ CL + cg-ca ^ CD + cg-caHT ^ CL_HT + cg-caHT ^ CD_HT = 0

        Si descompones este momento, habra elementos que:

           Dependan de alfa: Formaran Cm_alpha
           Dependan de nada: formaran Cm0
           Dependan de delta_e: Lo incluimos
           Dependeran de q:  Lo incluimos?
           Dependeran de d(alfa): Lo incluimos?  Con teoria del retardo




        cg-ca = SM

        cg-caHT = g.lv

        za =

        L_HT =    ( -0.0068 + 0.7798 *   ALPHA_HT   )  * (PRESION DINAMICA COLA)  /  (PRESION DINAMICA AVION)

        D_HT =    ( (POLAR) * ALPHA _HT  )   * (PRESION DINAMICA COLA)  /  (PRESION DINAMICA AVION)
    """


    """
    ATR
    xcg: 12.41 measured from tip of aircraft, in OpenVSP. You need to change in the files
        Thesis Hamburgo says Xcg = 11.5586, and allowable Xcg between 0.14 and 0.27 CMAC 
        ATR72-500 manual says Xcg = (11.4743-12.096). Maybe we should change OpenVSP value. And also forward the wing
        to 11.242
            c.dg. = 11.75
            c.dg.s without horizontal tail 11.571
            zcdg = 1.315
            zdg without horizontal tail 1.262
         
    Xlew: distance from tip of aircraft to leading edge of the wing in the root. In OpenVSP is 11.38 m.
          In ATR72-600 measured manually, this distance is of 11.242 meters. It appears MAC match with the
          chord in the wing root. 
    
    CMAC: Mean aerodynamic chord length. ATR72-500 manual says 2.303 . OpenVSP says 2.324481
    
    XLEMAC: Distance from tip of aircraft to leading edge of MAC. In ATR72-500 (this is an ATR72-600)
           is 11.24. So the center of gravity is around 10 - 39% of MAC, this is, from the tip, 11.4743 to 12.14217 m 
           
           Assuming that MAC is situated so that its 25% point is the MAC of the aircraft, then MAC is 11.82. 
           In OPENVSP however, the neutral point is calculated in x=12.89 
           If you just run the wing, then it is given around 12.07.
    
    
    
    NP: Neutral Point.  For the wing leading edge in the root at 11.38 m OpenVSP calculates it in 12.8213492 m 
       Its more backwards as there is an horizontal stabilizer.
       
       
    
        
        
       ATR
    
          5.2458*((11.75-11.23)/2.3244810)-(0.780797 * 60.8477 / 11.13)*(1-0.247779354)*1.0435*((11.13*14.09)/(60.8477*2.3244810)) = -2.5415
          
          El valor dado por OpenVSP es 2.47605 asi que bien

  
          OPENVSP DA -2.69  LUEGO MUY MUY BIEN !!!
    
          Cl_alpha_wingbody = 5.2458
          x_cg = 11.75
          x_ca_wb = 11.23
          CMA = 2.3244
          
          nt * (1 - d_eps/d_alpha) = es el necesario para que 0.7808 sea 0.7141.  Como d_eps / d_alpha es 0.2477 podemos saber que nt = 1.0435
          
          Vt (coef volumen) = Sh * lh  /  Sw * CMA   Sw = 60.8477  Sh = 11.13   lh = 25.84 - 11.75 = 14.09
          
          at = 0.7808 * Sw / Sh   =  4.1046
          
       


    
    
    In Development of a new methodology for the prediction of aircraft fuselage aerodynamic characteristics, Vincenzo
    Cusati, he calculated, just for the fuselage of ATR72, several coefficients, you can have it also for other things
    
    CD0_fuselage = 0.0085
    CM (α=0) = −0.0832
    CMα = 0.0222
     
    """


    """
    DECOL
    
    xcg:   x = 713 mm from tip of aircraft
           z = 0 mm from OpenVSP origin. 75 mm below propellers.
          
           Without horizontal tail --> x = 686 mm from tip of aircraft
                                       z = -2 mm
          
       
    Xlew: distance from tip of aircraft to leading edge of the wing in the root. In OpenVSP is 610.124 mm.
    
    
    CMAC: Mean aerodynamic chord length. OpenVSP says 250mm
    
    XLEMAC: Distance from tip of aircraft to leading edge of MAC. In OpenVSP 610.124 mm
           Center of gravity is 41.15% of MAC, this is, from the leading edge 103 mm, from the tip 713 mm.
     
    
    NP: Neutral Point. Given Static Margin (Pg 274 Eric's thesis) : 15%; Distance from leading edge = 140.5. Distance from tip 750.5
        OpenVSP gives 762.6485872 from tip SM = 0.1985943
    
        
        
       ATR
    
          4.4587*((0.713-0.64498)/0.25)-(0.7049 * 0.5 / 0.0877)*(1-0.2891)*1.03908854822*((0.0877*1.096)/(0.5*0.25)) = -1.06963098975


          OPENVSP DA -1.001088  LUEGO MUY MUY BIEN !!!
    
          Cl_alpha_wingbody = 4.4587
          x_cg = 0.713
          x_ca_wb = 0.64498
          CMA = 0.25
          
          nt * (1 - d_eps/d_alpha) = es el necesario para que 0.7049 sea 0.5267.  Como d_eps / d_alpha es 0.28091 podemos saber que nt = 1.03908854822

          Vt (coef volumen) = Sh * lh  /  Sw * CMA   Sw = 0.5  Sh =0.0877    lh = 1.809 - 0.713 = 1.096
          
          at = 0.7049 * Sw / Sh  
          
       
    c.dg. = 
    c.dg.s without horizontal tail
    zcdg = 
    zdg without horizontal tail 
   
    

    """

    return Cm, CL_tail


