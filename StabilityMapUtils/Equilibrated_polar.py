"""
File for calculating an equilibrated polar.
We are in longitudinal and in equilibrium.

Variables are : H, V, alpha, theta, gamma, delta_e, delta_x
Equations are: 2 forces, 1 pitch moment, 1 theta = alpha + gamma

You have to fix 3 variables to find solution.

    1 varible to fix: H.

    2 variable to fix: V

    3 variables to fix: alpha (-5, 20) every degree


Calculating CL and Cd for equilibrium

Vamos a hacer por el momento que te de la polar dandole una H y una velocidad y una horquilla de alphas

"""


import numpy as np
import math
from scipy.optimize import minimize
from scipy.optimize import least_squares
import sys
from StabilityMapUtils import Longitudinal
from StabilityMapUtils import AeroForces
import matplotlib.pyplot as plt




def Equilpolar(Matrix, CoefMatrix, Mach, g, PW):

    H = 762
    V = 40.145
    step = 1
    #flaps??!

    atmo = g.GetAtmo(H)
    M_base = V/atmo[0]

    Coef_base = AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)
    g.Matrix_no_tail_terms = AeroForces.CoefInterpol(M_base, Matrix[:, 1:], Mach)
    # Recalculate PW: No need; it does not depend on V nor on H

    # Long_variables;   (to define a combination of 2 or less)

    alpha_long = np.linspace(-5*math.pi/180, 15*math.pi/180, num=27)
    #alpha_long = np.array([1*np.pi/180])
    gamma_long = "TS"
    V_long = V
    de_long = "TS"
    dx_long = "TS"
    theta_long = "TS"



    k_list = []
    x_list = []
    fixtest_list = []
    diccons_list = []
    CL = []
    CD = []
    Cm = []
    dx_array = np.zeros(len(alpha_long))
    de_array = np.zeros(len(alpha_long))


    for i in range(len(alpha_long)):

           vars = [alpha_long[i], gamma_long, V_long, de_long, dx_long, theta_long]

           k, x, fixtest, diccons = Longitudinal.Long_Equilibrium(Coef_base, atmo, g, PW, vars)

           k_list.append(k)
           x_list.append(x)
           fixtest_list.append(fixtest)
           diccons_list.append(diccons)

           V_vect = np.ones(g.N_eng) * V
           Fx_vec = g.Thrust(np.full(g.N_eng, x[-1]), V_vect, atmo)
           Tc = Fx_vec/(2*atmo[1]*g.Sp*V**2)

           F = AeroForces.CalcForce_aeroframe_DEP(V, Coef_base, np.array([x[0], 0, 0, 0, 0, 0, x[7], 0]), Tc, atmo, g, PW)

           CL.append((-F[2])/(0.5*atmo[1]*V**2*g.S))
           CD.append((-F[0])/(0.5*atmo[1]*V**2*g.S))
           Cm.append((F[4])/(0.5*atmo[1]*V**2*g.S*g.c))

           de_array[i] = x[7]
           dx_array[i] = x[-1]


    CL_array = np.array(CL)
    CD_array = np.array(CD)
    Cm_array = np.array(Cm)



    fig1 = plt.figure()
    ax1 = fig1.gca()
    ax1.plot(alpha_long*180/np.pi, CL_array, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
    ax1.set_xlabel('alpha (°)')
    ax1.set_ylabel('CL')
    ax1.legend()
    ax1.grid()
    fig1.tight_layout()


    fig2 = plt.figure()
    ax2 = fig2.gca()
    ax2.plot(alpha_long*180/np.pi, CD_array, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
    ax2.set_xlabel('alpha (°)')
    ax2.set_ylabel('CD')
    ax2.legend()
    ax2.grid()
    fig2.tight_layout()


    fig3 = plt.figure()
    ax3 = fig3.gca()
    ax3.plot(CD_array, CL_array, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
    ax3.set_xlabel('CD')
    ax3.set_ylabel('CL')
    ax3.legend()
    ax3.grid()
    fig3.tight_layout()


    fig4 = plt.figure()
    ax4 = fig4.gca()
    ax4.plot(alpha_long*180/np.pi, Cm_array, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
    ax4.set_xlabel('alpha (°)')
    ax4.set_ylabel('Cm')
    ax4.legend()
    ax4.grid()
    fig4.tight_layout()

    fig5 = plt.figure()
    ax5 = fig5.gca()
    ax5.plot(alpha_long*180/np.pi, dx_array, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
    ax5.set_xlabel('alpha (°)')
    ax5.set_ylabel('dx')
    ax5.legend()
    ax5.grid()
    fig5.tight_layout()


    fig6 = plt.figure()
    ax6 = fig6.gca()
    ax6.plot(alpha_long*180/np.pi, de_array*180/np.pi, label="$α_i$, $T_c$ = {0:0.3f}".format(Tc[0]), linestyle='-.', color='0.25')
    ax6.set_xlabel('alpha (°)')
    ax6.set_ylabel('de')
    ax6.legend()
    ax6.grid()
    fig6.tight_layout()

    plt.show(block=True)   # added to plot correctly




    return k_list, x_list, fixtest_list, diccons_list



