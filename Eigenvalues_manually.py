
"""

author: david.planas-andres

Module extracts information from eigenvalues, for instance, not damped natural frequency, damping coefficient or time constant.

Module calculates "analytically" the dimensional eigenvalues using the aerodynamic coefficients


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



def Eig_info(Longieigvals,Lateigvals,g,fixtest):

       V=fixtest[0]

       eigenvalues = np.concatenate((Longieigvals,Lateigvals))
       matrix = np.zeros(   (  (len(eigenvalues)), (3)  )    )
       matrix_adim = np.zeros(   (  (len(eigenvalues)), (3)  )    )


       for i in range(len(eigenvalues)):

            if abs(eigenvalues[i].imag) >= 10e-9:

               matrix[i,0]=((eigenvalues.real[i])**2+(eigenvalues.imag[i])**2 )**0.5
               matrix[i,1]=(-eigenvalues.real[i])/ (((eigenvalues.real[i])**2+(eigenvalues.imag[i])**2 )**0.5)

               matrix_adim[i,1]=matrix[i,1]
               if i < 4:
                  matrix_adim[i,0]=matrix[i,0]*(g.c/(2*V))
               else:
                  matrix_adim[i,0]=matrix[i,0]*(g.b/(2*V))



            else:

               matrix[i,2]=-1/eigenvalues.real[i]

               if i <4:
                  matrix_adim[i,2]=-1/eigenvalues.real[i]*(g.c/(2*V))
               else:
                  matrix_adim[i,2]= -1/eigenvalues.real[i]*(g.b/(2*V))




       return matrix, matrix_adim


# rows: longitudinal eigenvalues, lateral eigenvalues
# colons: natural frequency, damping coefficient, time constant.
# adim stands for adimensionalized














def Eigenvalues(Aero_Derivatives_adim,x,fixtest,atmo,g,PW,CoefMatrix):


          rho = atmo[1]

          V = fixtest[0]
          alpha = x[0]
          beta = fixtest[1]
          gamma = fixtest[2]
          omega = fixtest[3]
          p = x[1]
          q = x[2]
          r = x[3]
          phi = x[4]
          theta = x[5]
          I = np.array([[g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz]])

          #for CLs and CDs calculus
          V_vect = np.ones(g.N_eng) * V * np.cos((-np.sign(g.PosiEng)) * beta + g.wingsweep) - r * g.PosiEng
          sub_vect = np.array([alpha, beta, p, q, r])
          if g.nofin == False:
              sub_vect = np.append(sub_vect, [x[6], x[7], x[8]])
          else:
              sub_vect = np.append(sub_vect, [x[6], x[7]])
          Fx_vec = g.Thrust(x[-g.N_eng:], V_vect)
          Tc = Fx_vec / (2 * rho * g.Sp * V_vect ** 2)
          F = AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), Tc, atmo, g, PW)



          mulong = g.m/(rho*g.S*g.c*0.5)
          mulat =g.m/(rho*g.S*g.b*0.5)




          C_Yr = Aero_Derivatives_adim[3, 1]
          C_Ybeta = Aero_Derivatives_adim[-3,1]

          C_Zs = - g.m * 9.81 * np.cos(theta) / (0.5*rho*V**2*g.S)
          C_Zu = Aero_Derivatives_adim[-4,2]
          C_Zq = Aero_Derivatives_adim[2,2]
          C_Zalfa = Aero_Derivatives_adim[0,2]
          C_Zdevalfa = 0


          C_mq = Aero_Derivatives_adim[2,4]
          C_malfa = Aero_Derivatives_adim[0,4]
          C_mdevalfa = 0




          C_lp = Aero_Derivatives_adim[1,3]
          C_lbeta = Aero_Derivatives_adim[-3,3]
          C_lr = Aero_Derivatives_adim[3, 3]


          C_nr = Aero_Derivatives_adim[3,5]
          C_nbeta = Aero_Derivatives_adim[-3,5]
          C_np = Aero_Derivatives_adim[1,5]


          I_x = I[0,0] / (rho*g.S*(g.b*0.5)**3)
          I_y = I[1,1] / (rho*g.S*(g.c*0.5)**3)
          I_z = I[2,2] / (rho*g.S*(g.b*0.5)**3)
          J_xz = -I[0,2] / (rho*g.S*(g.b*0.5)**3)








          C_Ls = F[2]/(0.5*rho*V**2*g.S)
          C_Ds = F[0]/(0.5*rho*V**2*g.S)





          #Phugoid:

          Wn_ph=((C_Zs*(2*C_Zs+C_Zu)/(2*mulong*(2*mulong+C_Zq)))**0.5) *(2*V/g.c)



          xi_ph=(2-(-1))/(2*(2**0.5)*(np.tan(alpha)+C_Ls/C_Ds))    #formula from isae

          lambda_ph= np.complex(-xi_ph*Wn_ph,Wn_ph*(1-xi_ph**2)**0.5)



          #Short_Period

          Wn_sp= (    (       (C_Zalfa*C_mq - (2*mulong+C_Zq)*C_malfa)  /  (I_y*(2*mulong-C_Zdevalfa))          )**0.5 )* (2*V/(g.c))



          xi_sp= -((2*mulong-C_Zdevalfa)*C_mq+C_Zalfa*I_y+(2*mulong + C_Zq)*C_mdevalfa)/( 2* ((2*mulong-C_Zdevalfa)*I_y*(C_Zalfa*C_mq - (2*mulong+C_Zq)*C_malfa) )**0.5)

          lambda_sp = np.complex(-xi_sp*Wn_sp,Wn_sp*(1-xi_sp**2)**0.5)


          #Roll subsidence

          lambda_rs= (C_lp/I_x) * (V/(g.b*0.5))



          #Spiral

          lambda_spiral=   (      (-C_Zs*(C_lbeta*C_nr-C_nbeta*C_lr))/   (    (2*mulat - C_Yr)*(C_nbeta*C_lp - C_lbeta*C_np) + C_Ybeta*(C_lp*C_nr-C_np*C_lr)     )  )* (V/(g.b*0.5))



          #Dutch Roll

          coef=[I_x*I_z,-(I_x*C_nr+I_z*C_lp),(I_x*C_nbeta+C_nr*C_lp-C_lr*C_np+J_xz*C_lbeta),-(C_nbeta*C_lp-C_lbeta*C_np)]

          roots = np.roots(coef)

          lambda_dr=roots[1] * (V/(g.b*0.5))



#          return np.hstack((Wn_ph,xi_ph,Wn_sp,xi_sp,lambda_rs,lambda_spiral,roots))

          return np.hstack((lambda_ph,lambda_sp,lambda_rs,lambda_spiral,lambda_dr))