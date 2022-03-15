
"""

author: david.planas-andres

Generation of a simple sample for the use of 1 Dimension splines for calculus of simple aerodynamic coefficients.


Vector x and fix to be assembled and varied
#    x = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i(hay 12), V , beta , gamma, omega] = x + fix

Function Constraints_DEP here gives the vector   CD Cy CL Cl Cm Cn

    COEF_MATRIX_DEF is a matrix of 3 dimension with dimension:  (number_points) rows x 6 colons x 24 variables (25 if rudder):
    (number_points) roxs for each of (number_points) interpolations of CL CD CM Cl Cn Cy inside each variable
    6 colons for CD Cy CL Cl Cm Cn
    24 times the 2D matrix for each variable of the "state vector" (25 if rudder)

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








def aero_coefficients(x, fix, CoefMatrix, atmo, g, PropWing):


     phimax = 5  # in degree the max bank angle authorized
     alphamax = 25  # in degree, stall bound for trimming
     deltaRmax = 30  # in degree
     ThrottleMax = 1  # max thrust level
     ThrottleMin = 1e-9  # min thruttle, don't accept 0 thrust
     number_points = 10





     #               alfa                            p            q           r                     phi                                     theta                        delta_a                                delta_e                                delta_r
     # bnds=( (-5*math.pi/180,alphamax*math.pi/180), (-0.2,0.2), (-0.2,0.2), (-0.2,0.2), (-phimax/180*math.pi,phimax/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-30/180*math.pi,30/180*math.pi), (-20/180*math.pi,20/180*math.pi), (-deltaRmax/180*math.pi,deltaRmax/180*math.pi))

     bnds2 = (( x[0] - 3 * math.pi / 180, x[0] + 3* math.pi / 180), (-0.2, 0.2), (-0.2, 0.2), (-0.2, 0.2),
           (-phimax / 180 * math.pi, phimax / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi),
           (-10 / 180 * math.pi, 10 / 180 * math.pi), ( x[7] -5 / 180 * math.pi, x[7] + 5 / 180 * math.pi),
           (-10 / 180 * math.pi, 10 / 180 * math.pi))




     if g.hangar['aircraft']=='ATR72':

        limfix = ( (fix[0]-0.9, fix[0]+0.9), (-8 / 180 * math.pi, 8 / 180 * math.pi), (-5 / 180 * math.pi, 5 / 180 * math.pi), (-0.2, 0.2))

     elif g.hangar['aircraft']=='DECOL':

        limfix=( (21,26),(-5/180*math.pi,5/180*math.pi),(-5/180*math.pi,5/180*math.pi),(-0.2,0.2) )







     bnds_eng=((ThrottleMin,ThrottleMax), (ThrottleMin,ThrottleMax))

     for i in range(int(g.N_eng/2)):
         bnds2=bnds2+bnds_eng

     bnds2=bnds2+limfix




     x = np.concatenate((x, fix))




     COEF_MATRIX = np.zeros( (number_points,6,len(x)) )
     VAR_MATRIX=np.zeros((number_points,len(x)))








     for i in range(len(x)):

         x_dev = list(x)
         x_dev = np.asarray(x_dev)

         for j in range(number_points):

              x_dev[i] = bnds2[i][0]+j*(bnds2[i][1]-bnds2[i][0])/number_points



              VAR_MATRIX[j, i] = x_dev[i]

              COEF_MATRIX[j, :, i] = Constraints_DEP(x_dev, CoefMatrix, atmo, g, PropWing)






     Int_list_1=[]
     Int_list_3=[]
     Coef_Dev_Interpolation = []
     Coef_Dev_Interpolation_2 = []


     Int_list_5=[]








     for i in range(len(x)):

          for j in range(6):

              Int_list_1.append(interp1d(VAR_MATRIX[:,i], COEF_MATRIX[:, j, i]))

              u=IUS(VAR_MATRIX[:, i], COEF_MATRIX[:, j, i])

              Int_list_3.append(u.derivative())
              



          Int_list_2=list(Int_list_1)
          Int_list_4=list(Int_list_3)
          Coef_Dev_Interpolation.append(Int_list_2)
          Coef_Dev_Interpolation_2.append(Int_list_4)
          Int_list_1.clear()
          Int_list_3.clear()





     plt.plot(VAR_MATRIX[:,0], COEF_MATRIX[:,2,0])



     Aero_Derivatives=np.zeros((len(x),6))

     for i in range(len(x)):
         for j in range(6):

             Aero_Derivatives[i,j]=Coef_Dev_Interpolation_2[i][j](x[i])





     #For adimensionalizing derivatives with respect to p,q,r,V

     b1 = [[Aero_Derivatives[x][y] for y in range(len(Aero_Derivatives[0]))] for x in range(len(Aero_Derivatives))]
     Aero_Derivatives_adim=np.asarray(b1)

     Aero_Derivatives_adim[1,:] = Aero_Derivatives_adim[1,:]*(2*fix[0]/g.b)
     Aero_Derivatives_adim[2,:] = Aero_Derivatives_adim[2,:]*(2*fix[0]/g.c)
     Aero_Derivatives_adim[3,:] = Aero_Derivatives_adim[3,:]*(2*fix[0]/g.b)
     Aero_Derivatives_adim[-4,:] = Aero_Derivatives_adim[-4,:] * (fix[0])







     return Aero_Derivatives, Aero_Derivatives_adim



"""

     VAR_MATRIX: A matrix with 25 colons (for each variable) and rows equal to number_points, for each variation in 
     the variables


     Coef_Dev_Interpolation: Is a list with 24 lists inside (una for each variable) Each list has inside as well 6 
     functions of interpolation (CL CD CM Cl Cn Cy)
     So the element Coef_Dev_Interpolation[2][4] is the interpolation function that corresponds to the calculation of
     Cn (element 4) when varying pitch rate q (element 2 of assembled vector x).
     Value of the approach of the element [i][j] in the point x can be calculated as  Coef_Dev_Interpolation[i][j](x)


     #Coef_Dev_Interpolation_2 is exactly the same but contains the derivatives of those interpolation functions
     #Derivative in element [i][j] in point x is calculated as  Coef_Dev_Interpolation_2[i][j](x)

"""











def Constraints_DEP(x, CoefMatrix, atmo, g, PropWing):


    #x = [alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i(hay 12), V , beta , gamma, omega] = x + fix

    rho = atmo[1]

    n_eng = int(g.N_eng / 2)

    # --- Now prepare variables for equations ---
    V = x[-4]
    alpha = x[0]
    beta = x[-3]
    gamma = x[-2]
    omega = x[-1]
    p = x[1]
    q = x[2]
    r = x[3]
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

















