"""

author: david.planas-andres

Systematic sample generation for orthogonal -least squares algorithm APRICOT

The following variables are varied for longitudinal and lateral cases
   x + fix = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_xi(hay 12), V , beta , gamma, omega]


 Given example. For a function y=f(x1,x2,x3) . In order to get a systematic sample working in apricott a 3 points
 sample is taken (0,1,2) for the three variables, so number of points = variations^variables = (3^3=27)
 Apricott requires the order of the sample to be as following (systematic)



     # x1: 0   0   0       0   0   0       O   O   O       1   1   1       1   1   1       1   1   1       2   2   2       2   2   2       2   2   2
     # x2: 0   0   0       1   1   1       2   2   2       0   0   0       1   1   1       2   2   2       0   0   0       1   1   1       2   2   2
     # x3: 0   1   2       0   1   2       0   1   2       0   1   2       0   1   2       0   1   2       0   1   2       0   1   2       0   1   2

     #It is required to build a matrix of:
     # rows: 3 (number of variables)
     # colons: 27 (variations^variables, 3^3)

     #VARIABLE 1: suffers a variation each  number of variations^2
     #VARIABLE 2: suffers a variation each  number of variations^1
     #VARIABLE 3: suffers a variation each  number of variations^0



     variables=3
     variations=3
     samplevector=np.zeros((variables,variations**variables))

     for i in range(variations):              #from 0 to 2
         for j in range(variations):          #from 0 to 2
               for k in range(variations):    #from 0 to 2

                    samplevector[variables-3,i*variations**2+j*variations+k] = i
                    samplevector[variables-2,i*variations**2+j*variations+k] = j
                    samplevector[variables-1,i*variations**2+j*variations+k] = k
"""


from StabilityMapUtils import AeroForces
import numpy as np
import math
import scipy.linalg
import scipy.io #input/output with matlab



def Sample_generation(x, fix, CoefMatrix, atmo, g, PropWing):





    if g.hangar['aircraft']=='ATR72':

        phimax = 5  # in degree the max bank angle authorized
        alphamax = 25  # in degree, stall bound for trimming
        deltaRmax = 30  # in degree
        ThrottleMax = 1  # max thrust level
        ThrottleMin = 1e-9  # min throttle, don't accept 0 thrust
        V = fix[0]

        #               alfa                            p            q           r                     phi                                     theta                        delta_a                                delta_e                                delta_r
        # bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-23/180*math.pi,13/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))



        bnds = (( x[0] - 3 * math.pi / 180, x[0] + 3* math.pi / 180), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))), (-0.2*(g.c/(2*V)), 0.2*(g.c/(2*V))), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))),
                (-phimax / 180 * math.pi, phimax / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi),
                (-10 / 180 * math.pi, 10 / 180 * math.pi), ( x[7] -5 / 180 * math.pi, x[7] + 5 / 180 * math.pi),
                (-10 / 180 * math.pi, 10 / 180 * math.pi))


        limfix = ( ((V-10)/V, (V+10)/V), (-5 / 180 * math.pi, 5 / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi), (-0.2, 0.2))








    elif g.hangar['aircraft']=='DECOL':

        phimax = 10  # in degree the max bank angle authorized
        alphamax = 25  # in degree, stall bound for trimming
        deltaRmax = 30  # in degree
        ThrottleMax = 1  # max thrust level
        ThrottleMin = 0.0001  # min throttle, don't accept 0 thrust
        V = fix[0]

        #               alfa                            p             q            r                   phi                                     theta                        delta_a                                delta_e                                delta_r
        #bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))




        bnds = ((0 * math.pi / 180, 5 * math.pi / 180), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))), (-0.2*(g.c/(2*V)), 0.2*(g.c/(2*V))), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))),
                (-phimax / 180 * math.pi, phimax / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi),
                (-10 / 180 * math.pi, 10 / 180 * math.pi), (x[7] -5 / 180 * math.pi, x[7] + 5 / 180 * math.pi),
                (-10 / 180 * math.pi, 10 / 180 * math.pi))


        limfix=( (21/V,26/V), (-5/180*math.pi,5/180*math.pi), (-5/180*math.pi,5/180*math.pi), (-0.2,0.2) )



    #For longitudinal
    bnds_eng_long = ((ThrottleMin, ThrottleMax))
    bnds_long = bnds+(bnds_eng_long,) + limfix



    x = np.concatenate((x, fix))
    # Adimensionalizing variables:

    x[1] = x[1] / (2 * V / g.b)         #p
    x[2] = x[2] / (2 * V / g.c)         #q
    x[3] = x[3] / (2 * V / g.b)         #r
    x[-4] = x[-4] / V                   #V






    #LONGITUDINAL
    # (CD,CL,Cm)
    #Variables (2): alpha ,delta_xi (all engines)

    #variables and position in vector x (starting on 0)
    # alpha   0
    # delta_xi: -(g.N_eng+4):-4


    #IN LONGITUDINAL, ALL ENGINES ARE VARIED IN THE SAME WAY, dx FROM 0 TO 1

    variations = 10
    variables = 2
    testvector = np.zeros((len(x), (variations**variables)))
    Xsample_longitudinal = np.zeros((variables, variations ** variables))


    Coefs=np.zeros((6, (variations**variables)))

    for i in range(variations):
        for j in range(variations):

                        testvector[:, i*variations+j] = x

                        testvector[0, i*variations+j] = bnds_long[0][0] +i*(bnds_long[0][1] - bnds_long[0][0])/ (variations-1)
                        #ENGINES
                        testvector[-(g.N_eng+4):-4, i*variations+j] = bnds_long[-5][0]+j*(bnds_long[-5][1] - bnds_long[-5][0])/(variations-1)


                        Coefs[:, i*variations+j] = Constraints_DEP(testvector[:, i*variations+j], CoefMatrix, atmo, g, PropWing,V)

    Xsample_longitudinal[0, :] = testvector[0, :]
    Xsample_longitudinal[1, :] = testvector[12, :]


    CD_sample = Coefs[0, :]
    CL_sample = Coefs[2, :]
    Cm_sample = Coefs[4, :]



    Xsample_longitudinal, CD_sample, CL_sample,Cm_sample = Python_to_Matlab(Xsample_longitudinal, CD_sample,CL_sample, Cm_sample)


    return Xsample_longitudinal, CD_sample, CL_sample, Cm_sample






def Python_to_Matlab(x_long, CD, CL,Cm):

    scipy.io.savemat('alpha_dx_figure.mat', dict(x_long=x_long, CD=CD, CL=CL, Cm=Cm))

    return x_long, CD, CL, Cm












def Constraints_DEP(x, CoefMatrix, atmo, g, PropWing,Vfix):


    rho = atmo[1]

    n_eng = int(g.N_eng / 2)

    # --- Now prepare variables for equations ---
    V = x[-4]*Vfix
    alpha = x[0]
    beta = x[-3]
    gamma = x[-2]
    omega = x[-1]

    p = x[1]*(2*Vfix/g.b)
    q = x[2]*(2*Vfix/g.c)
    r = x[3]*(2*Vfix/g.b)

    phi = x[4]
    theta = x[5]
    I = np.array([[g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz]])



    # --- Compute aerodynamic forces ---
    # here subvector  must be : (alpha, beta, p, q, r, da, de,dr,  dx)
    sub_vect = np.array([alpha, beta, p, q, r])
    if g.nofin == False:
        sub_vect = np.append(sub_vect, [x[6], x[7], x[8]])  # rudder is allowed
    else:
        sub_vect = np.append(sub_vect, [x[6], x[7]])  # no fin allowed, default case


    V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng
    Fx_vec = g.Thrust(x[-(g.N_eng+4):-4], V_vect)
    Tc = Fx_vec / (2 * rho * g.Sp * V_vect ** 2)

    F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PropWing)
    # F contains forces and moments in wind reference system, just the aerodynammic, not thrust

    F[0] = np.abs(F[0])

    F[2] = np.abs(F[2])

    Coefs=np.zeros(len(F))

    for i in range(len(F)):
        if i == 0 or i == 1 or i == 2:
            Coefs[i] = F[i] / (0.5 * rho * V**2 * g.S)

        elif i == 4:
            Coefs[i] = F[i] / (0.5 * rho * V**2 * g.S * g.c)

        else:
            Coefs[i] = F[i] / (0.5 * rho * V ** 2 * g.S * g.b)


    return Coefs