"""
Using Treffz plane to calculate the induced drag
https://www.youtube.com/watch?v=zmz2gV3HHpU&ab_channel=BYUFLOWLab

We need to retrieve the z position of each panel. Better if we do this at the some time than the other magnitudes.
Opening and reading files is very high time consuming.
"""


import numpy as np

def Trefftz_drag(SortedCoef, plane,V):
    """
    There are N panels in the wing. A magnitude can be computed in the middle point of each panel,
    middle points: i = 1,...,N   (python i = 0,...,N-1)

    There are N+1 nodal points at the edges of each panel.A magnitude can also be computed in these
    nodal points: j = 1,..., N+1   (python i = 0,...,N)

    Variables       Points of computation
    -------------------------------------
    GAMMA               middle i
    gamma               nodal j
    V                   middle i
    y                   nodal j
    y_hat               middle i
    z                   nodal j
    z_hat               middle i

    For understanding the development below is is strongly recommended to watch:
    https://www.youtube.com/watch?v=zmz2gV3HHpU&list=LLCgI17svp5Ma_uOz3CHUvdA&index=4&ab_channel=BYUFLOWLab


    Induced drag is:     Di = 0.5 * rho * sum(Vn * GAMMA * S) , magnitudes computed at the center of panels.

    Gamma_i = 0.5 * chord * V * Cl


            = -Gamma_1                  if j=1
    gamma_j =  Gamma_(i-1) - Gamma_i    if j=2,...,n
            =  Gamma_n                  if j=n+1


    Inputs: SortedCoef
            plane
    """


    GAMMA = 0.5 * SortedCoef['LocalChord'] * SortedCoef['Cl'] * SortedCoef['Vep_total']

    gamma = np.hstack((-GAMMA[0], -np.diff(GAMMA), GAMMA[-1]))


    y_hat = SortedCoef['Yposi']
    y = DiffclPosi = np.hstack(((-plane.b/2), SortedCoef['Yposi'][1:] - np.diff(SortedCoef['Yposi'])/2, (plane.b/2)))

    z_hat = np.zeros(len(y_hat))
    z = np.zeros(len(y))

    dY = abs(y_hat[1]-y_hat[2])


    Di_i = np.zeros(len(GAMMA))
    for i in range(len(GAMMA)):

        Di_i[i] = (1/4*np.pi) * np.sum(GAMMA[i]*gamma * ((y-y_hat[i])*(y[i+1]-y[i]) + (z-z_hat[i])*(z[i+1]-z[i]))/(((y-y_hat[i])**2 + (z-z_hat[i])**2)))

    Di = np.sum(Di_i * dY)/(0.5*V**2*plane.S)



    """
    This did not work since the sharps distributions of lift due to "artificial" augmentation of lift due to blowing 
    and/or flaps lead to not realistic circulations in some point, so gamma  = Gamma1 - Gamma2 is very unrealistic and
    usually the slope is REALLY well captured but there is a vertical desplacement of the curve.
    
    For induced drag: 
       U can implement Jamesson
       U can leave LLT drag minus a constant (lets say, with a correction)
    
    Independently u can try to improve the post-stall, but post stall allows you not to be able to trim the aircraft in stall
    so is not hat bad. 
    """


    return Di




