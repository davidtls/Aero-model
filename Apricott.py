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

import numpy as np
import math
from StabilityMapUtils import AeroForces
from numpy.linalg import inv
import ReadFileUtils as Read  # utils to read Xfoil file
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
import matplotlib.pyplot as plt


def Sample_generation(x, fix, CoefMatrix, atmo, g, PropWing):





     if g.hangar['aircraft']=='ATR72':

        phimax = 5  # in degree the max bank angle authorized
        alphamax = 25  # in degree, stall bound for trimming
        deltaRmax = 30  # in degree
        ThrottleMax = 1  # max thrust level
        ThrottleMin = 1e-9  # min thruttle, don't accept 0 thrust
        V=fix[0]

        #               alfa                            p            q           r                     phi                                     theta                        delta_a                                delta_e                                delta_r
        # bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))



        bnds = (( x[0] - 3 * math.pi / 180, x[0] + 3* math.pi / 180), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))), (-0.2*(g.c/(2*V)), 0.2*(g.c/(2*V))), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))),
             (-phimax / 180 * math.pi, phimax / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi),
             (-10 / 180 * math.pi, 10 / 180 * math.pi), ( x[7] -5 / 180 * math.pi, x[7] + 5 / 180 * math.pi),
             (-10 / 180 * math.pi, 10 / 180 * math.pi))


        limfix = ( ((V-0.9)/V, (V+0.9)/V), (-5 / 180 * math.pi, 5 / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi), (-0.2, 0.2))








     elif g.hangar['aircraft']=='DECOL':

        phimax = 10  # in degree the max bank angle authorized
        alphamax = 25  # in degree, stall bound for trimming
        deltaRmax = 30  # in degree
        ThrottleMax = 1  # max thrust level
        ThrottleMin = 0.0001  # min throttle, don't accept 0 thrust
        V=fix[0]

       #               alfa                            p             q            r                   phi                                     theta                        delta_a                                delta_e                                delta_r
       #bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))




        bnds = ((0 * math.pi / 180, 5 * math.pi / 180), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))), (-0.2*(g.c/(2*V)), 0.2*(g.c/(2*V))), (-0.2*(g.b/(2*V)), 0.2*(g.b/(2*V))),
             (-phimax / 180 * math.pi, phimax / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi),
             (-10 / 180 * math.pi, 10 / 180 * math.pi), ( x[7] -5 / 180 * math.pi, x[7] + 5 / 180 * math.pi),
             (-10 / 180 * math.pi, 10 / 180 * math.pi))


        limfix=( (21/V,26/V),(-5/180*math.pi,5/180*math.pi),(-5/180*math.pi,5/180*math.pi),(-0.2,0.2) )



     #For longitudinal
     bnds_eng_long=((ThrottleMin,ThrottleMax))
     bnds_long=bnds+(bnds_eng_long,) + limfix

     #For Lateral
     bnds_eng_lat = tuple()
     bnds_eng_lat=((x[10]-0.25,x[10]+0.25))
     bnds_lat=bnds+(bnds_eng_lat,) + limfix


     x = np.concatenate((x, fix))
     # Adimensionalizing variables:

     x[1] = x[1] / (2 * V / g.b)         #p
     x[2] = x[2] / (2 * V / g.c)         #q
     x[3] = x[3] / (2 * V / g.b)         #r
     x[-4] = x[-4] / V                   #V








     #LONGITUDINAL
     # (CD,CL,Cm)
     #Variables (5): alpha ,q , delta_e ,delta_xi (all engines) , V

     #variables and position in vector x (starting on 0)
     # alpha   0
     # q   2
     # delta_e   7
     # delta_xi: -(g.N_eng+4):-4
     # V  -4

     #IN LONGITUDINAL, ALL ENGINES ARE VARIED IN THE SAME WAY, dx FROM 0 TO 1


     variations=5
     variables=5
     testvector=np.zeros((len(x),(variations**variables)))
     Xsample_longitudinal = np.zeros((variables, variations ** variables))

     CD_sample = np.zeros((variations**variables))
     CL_sample = np.zeros((variations**variables))
     Cm_sample = np.zeros((variations**variables))

     Coefs=np.zeros((6,(variations**variables)))

     for i in range(variations):
         for j in range(variations):
               for k in range(variations):
                    for p in range(variations):
                        for q in range(variations):


                            testvector[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q]=x


                            testvector[0, i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds_long[0][0] +i*(bnds_long[0][1] - bnds_long[0][0])/ (variations-1)
                            testvector[2, i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds_long[2][0] +j*(bnds_long[2][1] - bnds_long[2][0])/(variations-1)
                            testvector[7, i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds_long[7][0] +k*(bnds_long[7][1] - bnds_long[7][0])/(variations-1)

                            #ENGINES
                            testvector[-(g.N_eng+4):-4, i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds_long[-5][0]+p*(bnds_long[-5][1] - bnds_long[-5][0])/(variations-1)

                            testvector[-4 ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds_long[-4][0]+q*(bnds_long[-4][1] - bnds_long[-4][0])/(variations-1)

                            Coefs[:, i*variations**4+j*variations**3+k*variations**2+p*variations+q] = Constraints_DEP(testvector[:, i*variations**4+j*variations**3+k*variations**2+p*variations+q], CoefMatrix, atmo, g, PropWing,V)

     Xsample_longitudinal[0, :] = testvector[0, :]
     Xsample_longitudinal[1, :] = testvector[2, :]
     Xsample_longitudinal[2, :] = testvector[7, :]
     Xsample_longitudinal[3, :] = testvector[12, :]
     Xsample_longitudinal[4, :] = testvector[-4, :]

     CD_sample = Coefs[0, :]
     CL_sample = Coefs[2, :]
     Cm_sample = Coefs[4, :]






























     #LATERAL
     # (Cy,Cl,Cn)
     # Variables (6): p  , r  , delta_a ,  delta_r  , delta_xi, beta

     # variables and position in vector x (starting on 0)
     # p   1
     # r   3
     # delta_a 6
     # delta_r 8
     # delta_xi: -(g.N_eng+4):-4
     # beta -3

     #IN LATERAL ENGINES SHOULD BE VARIED DIFFERENTLY IN EACH WING. WHILE ONES ARE GIVEN MORE THRUST,
     #THE OTHERS SHOULD BE GIVEN LESS THRUST
     # left wing (with y negative): x[(-g.N_eng-4):-(g.N_eng//2+4)] = [ -16  : -10 ]  (ATR)
     # right wing (with y positive): x[-(g.N_eng//2+4):-4]  =  [-10 : -4 ]  (ATR)
     #        bnds_lat=((2*x[9]-1,1), (1,2*x[9]-1))

     variations2 = 5
     variables2 = 6
     testvector2 = np.zeros((len(x), (variations2 ** variables2)))

     Xsample_lateral = np.zeros((variables2, variations2 ** variables2))

     CY_sample = np.zeros((variations2**variables2))
     Cl_sample = np.zeros((variations2**variables2))
     Cn_sample = np.zeros((variations2**variables2))

     Coefs2 = np.zeros((6, (variations2 ** variables2)))

     for i in range(variations2):
         for j in range(variations2):
               for k in range(variations2):
                    for p in range(variations2):
                        for q in range(variations2):
                            for r in range(variations2):


                                testvector2[:, i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r]=x


                                testvector2[1, i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r] = bnds_lat[1][0] +i*(bnds_lat[1][1]  - bnds_lat[1][0])/(variations2-1)
                                testvector2[3, i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r] = bnds_lat[3][0] +j*(bnds_lat[3][1] - bnds_lat[3][0])/(variations2-1)
                                testvector2[6, i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r] = bnds_lat[6][0] +k*(bnds_lat[6][1] - bnds_lat[6][0])/(variations2-1)
                                testvector2[8, i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r] = bnds_lat[8][0]+p*(bnds_lat[8][1] - bnds_lat[8][0])/(variations2-1)


                                if g.hangar['aircraft'] == 'ATR72':
                                    # ENGINES
                                    # Left
                                    testvector2[(-g.N_eng - 4):-(g.N_eng // 2 + 4), i * variations2 ** 5 + j * variations2 ** 4 + k * variations2 ** 3 + p * variations2 ** 2 + q * variations2 + r] = bnds_lat[-5][0] + q * (bnds_lat[-5][1] - bnds_lat[-5][0]) / (variations2 - 1)
                                    # Right
                                    testvector2[-(g.N_eng // 2 + 4):-4, i * variations2 ** 5 + j * variations2 ** 4 + k * variations2 ** 3 + p * variations2 ** 2 + q * variations2 + r] = bnds_lat[-5][1] + q * (bnds_lat[-5][0] - bnds_lat[-5][1]) / (variations2 - 1)

                                elif g.hangar['aircraft'] == 'DECOL':

                                    # ENGINES    EN DECOL SOLO VARIAN LOS DOS MOTORES EXTERIORES. EL RESTO TIENEN QUE TENER UN VALOR NOMINAL DE POTENCIA, CUAL? LES PONGO 0.5? los dejo a su valor de trimado
                                    # Left
                                    testvector2[(-g.N_eng - 4):(-g.N_eng -2), i * variations2 ** 5 + j * variations2 ** 4 + k * variations2 ** 3 + p * variations2 ** 2 + q * variations2 + r] = bnds_lat[-5][0] + q * (bnds_lat[-5][1] - bnds_lat[-5][0]) / (variations2 - 1)
                                    # Right
                                    testvector2[-6:-4, i * variations2 ** 5 + j * variations2 ** 4 + k * variations2 ** 3 + p * variations2 ** 2 + q * variations2 + r] = bnds_lat[-5][1] + q * (bnds_lat[-5][0] - bnds_lat[-5][1]) / (variations2 - 1)

                                testvector2[-3, i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r] = bnds_lat[-3][0]+r*(bnds_lat[-3][1] - bnds_lat[-3][0])/(variations2-1)



                                Coefs2[:,i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r]=Constraints_DEP(testvector2[:,i*variations2**5+j*variations2**4+k*variations2**3+p*variations2**2+q*variations2+r], CoefMatrix, atmo, g, PropWing,V)





     Xsample_lateral[0, :] = testvector2[1, :]
     Xsample_lateral[1, :] = testvector2[3, :]
     Xsample_lateral[2, :] = testvector2[6, :]
     Xsample_lateral[3, :] = testvector2[8, :]
     Xsample_lateral[4, :] = testvector2[10, :]
     Xsample_lateral[5, :] = testvector2[-3, :]


     CY_sample = Coefs2[1, :]
     Cl_sample = Coefs2[3, :]
     Cn_sample = Coefs2[5, :]





     """ 

     #CD
     #Variables: alpha, delta_e, deltax_i , V

     #variables and position in vector x (starting on 0)
     # alpha   : 0
     # delta_e : 7
     # deltax_i r  : 12
     # V       : -4

     variations=10
     variables=4
     testvector=np.zeros((len(x),(variations**variables)))
     CD_Xsample = np.zeros((variables, variations ** variables))
     CD_Ysample=np.zeros(variations**variables)
     Coefs=np.zeros((6,(variations**variables)))

     for i in range(variations):
         for j in range(variations):
               for k in range(variations):
                    for p in range(variations):

                       testvector[:,i*variations**3+j*variations**2+k*variations+p]=x

                       testvector[0  ,i*variations**3+j*variations**2+k*variations+p] = bnds2[0][0] +i*(bnds2[0][1]  - bnds2[0][0] )/variations
                       testvector[7  ,i*variations**3+j*variations**2+k*variations+p] = bnds2[7][0] +j*(bnds2[7][1]  - bnds2[7][0] )/variations
                       testvector[12 ,i*variations**3+j*variations**2+k*variations+p] = bnds2[12][0]+k*(bnds2[12][1] - bnds2[12][0])/variations
                       testvector[-4 ,i*variations**3+j*variations**2+k*variations+p] = bnds2[-4][0]+p*(bnds2[-4][1] - bnds2[-4][0])/variations

                       Coefs[:,i*variations**3+j*variations**2+k*variations+p]=Constraints_DEP(testvector[:,i*variations**3+j*variations**2+k*variations+p], CoefMatrix, atmo, g, PropWing)

     CD_Xsample[0,:]=testvector[0 ,:]
     CD_Xsample[1,:]=testvector[7 ,:]
     CD_Xsample[2,:]=testvector[12 ,:]
     CD_Xsample[3,:]=testvector[-4 ,:]
     CD_Ysample=Coefs[0,:]



     #CONSIDERACIONES
     #HAY QUE TENER CUIDADO CON TESTVECTOR, PORQUE SE LLAMA EN TODOS LOS BUCLES IGUAL, REVISA SI CUANDO DENTRO DEL BUCLE LE DICES QUE ES IGUAL QUE X SE REINICIA Y
     #TOMA EL VALOR DEL VECTOR DE ESTADO
     #PARA PONER TODOS LOS MOTORES DE UNA VEZ YO CREO QUE SERÍA:
     #      testvector[-(g.N_eng+4):-4 ,LO QUE SEA] = bnds2[-(g.N_eng+4):-4][0]+(LO QUE SEA)*(bnds2[12][1] - bnds2[12][0])/variations

     #Cy
     #Variables: p ,r , delta_r (if true) ,CT_i

     #variables and position in vector x (starting on 0)
     # p   1
     # r   3
     # delta r   8
     # CT_i
     # beta  -3


     variations=10
     variables=5
     testvector=np.zeros((len(x),(variations**variables)))
     CY_Xsample = np.zeros((variables, variations ** variables))
     CY_Ysample=np.zeros(variations**variables)
     Coefs=np.zeros((6,(variations**variables)))

     for i in range(variations):
         for j in range(variations):
               for k in range(variations):
                    for p in range(variations):
                        for q in range(variations):


                            testvector[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q]=x


                            testvector[1  ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[1][0] +i*(bnds2[1][1]  - bnds2[1][0] )/variations
                            testvector[3  ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[3][0] +j*(bnds2[3][1]  - bnds2[3][0] )/variations
                            testvector[8  ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[8][0] +k*(bnds2[8][1] - bnds2[8][0])/variations
                            testvector[12 ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[12][0]+p*(bnds2[12][1] - bnds2[12][0])/variations

                            testvector[-3 ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[-3][0]+q*(bnds2[-3][1] - bnds2[-3][0])/variations

                            Coefs[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q]=Constraints_DEP(testvector[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q], CoefMatrix, atmo, g, PropWing)

     CY_Xsample[0,:]=testvector[1 ,:]
     CY_Xsample[1,:]=testvector[3 ,:]
     CY_Xsample[2,:]=testvector[8 ,:]
     CY_Xsample[3,:]=testvector[12 ,:]
     CY_Xsample[4, :] = testvector[-3, :]
     CY_Ysample=Coefs[1,:]













     #CL
     #Variables  , alfa ,q , delta_e ,CT_i , V

     #variables and position in vector x (starting on 0)
     # alfa   0
     # q   2
     # delta_e   6
     # CT_i  12
     # V  -4


     variations=10
     variables=5
     testvector=np.zeros((len(x),(variations**variables)))
     Xsample = np.zeros((variables, variations ** variables))
     Ysample=np.zeros(variations**variables)
     Coefs=np.zeros((6,(variations**variables)))

     for i in range(variations):
         for j in range(variations):
               for k in range(variations):
                    for p in range(variations):
                        for q in range(variations):


                            testvector[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q]=x


                            testvector[0  ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[0][0] +i*(bnds2[0][1]  - bnds2[0][0] )/variations
                            testvector[2  ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[2][0] +j*(bnds2[2][1]  - bnds2[2][0] )/variations
                            testvector[6  ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[6][0] +k*(bnds2[6][1] - bnds2[6][0])/variations
                            testvector[12 ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[12][0]+p*(bnds2[12][1] - bnds2[12][0])/variations

                            testvector[-4 ,i*variations**4+j*variations**3+k*variations**2+p*variations+q] = bnds2[-4][0]+q*(bnds2[-4][1] - bnds2[-4][0])/variations

                            Coefs[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q]=Constraints_DEP(testvector[:,i*variations**4+j*variations**3+k*variations**2+p*variations+q], CoefMatrix, atmo, g, PropWing)

     Xsample[0,:] =testvector[0 ,:]
     Xsample[1,:] =testvector[2 ,:]
     Xsample[2,:] =testvector[6 ,:]
     Xsample[3,:] =testvector[12 ,:]
     Xsample[4,:] = testvector[-4, :]
     Ysample=Coefs[2,:]









     #Cl  y Cn match in variables, so we do just one time
     # Variables: p  , r  , delta_a ,  delta_r  , CT_i ,  beta

     # variables and position in vector x (starting on 0)
     # p   0
     # r   2
     # delta_a 6
     # delta_r 8
     # CT_i  12
     #beta -3

     variations=10
     variables=6
     testvector=np.zeros((len(x),(variations**variables)))
     Xsample = np.zeros((variables, variations ** variables))
     Ysample=np.zeros(variations**variables)
     Coefs=np.zeros((6,(variations**variables)))

     for i in range(variations):
         for j in range(variations):
               for k in range(variations):
                    for p in range(variations):
                        for q in range(variations):
                            for r in range(variations):


                                testvector[:,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r]=x


                                testvector[0  ,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r] = bnds2[0][0] +i*(bnds2[0][1]  - bnds2[0][0] )/variations
                                testvector[2  ,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r] = bnds2[2][0] +j*(bnds2[2][1]  - bnds2[2][0] )/variations
                                testvector[6  ,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r] = bnds2[6][0] +k*(bnds2[6][1] - bnds2[6][0])/variations
                                testvector[8  ,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r] = bnds2[8][0]+p*(bnds2[8][1] - bnds2[8][0])/variations

                                testvector[12 ,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r] = bnds2[12][0]+q*(bnds2[12][1] - bnds2[12][0])/variations
                                testvector[-3, i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r] = bnds2[-3][0]+q*(bnds2[-3][1] - bnds2[-3][0])/variations



                                Coefs[:,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r]=Constraints_DEP(testvector[:,i*variations**5+j*variations**4+k*variations**3+p*variations**2+q*variations+r], CoefMatrix, atmo, g, PropWing)

     Xsample[0,:] =testvector[0 ,:]
     Xsample[1,:] =testvector[2 ,:]
     Xsample[2,:] =testvector[6 ,:]
     Xsample[3,:] =testvector[8 ,:]
     Xsample[4,:] =testvector[12, :]
     Xsample[5,:] =testvector[-3,:]

     Ysample=Coefs[3,:]
     Ysample2 = Coefs[5, :]






     #Cm


"""






     return Xsample_longitudinal, Xsample_lateral, CD_sample, CY_sample, CL_sample, Cl_sample, Cm_sample, Cn_sample


















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


     Coefs=np.zeros(len(F))

     for i in range(len(F)):
        if i==0 or i==1 or i==2:
            Coefs[i] = F[i] / (0.5 * rho * V**2 * g.S)

        elif i==4:
            Coefs[i] = F[i] / (0.5 * rho * V**2 * g.S * g.c)

        else:
            Coefs[i] = F[i] / (0.5 * rho * V ** 2 * g.S *  g.b)


     return Coefs
























































