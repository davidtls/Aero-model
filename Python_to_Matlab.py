

"""

author: david.planas-andres

Stores the input sample and the output sample into a Matlab file to be used in Matlab


"""

import numpy as np
import math
import scipy.linalg
import scipy.io #input/output with matlab
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from scipy.optimize import  minimize
import sys




def Python_to_Matlab(x_long,x_lat, CD, CY, CL, Cl,Cm, Cn):

    scipy.io.savemat('test.mat', dict(x_long=x_long, x_lat=x_lat, CD=CD, CY=CY, CL=CL, Cl=Cl, Cm=Cm, Cn=Cn))

    return x_long,x_lat, CD, CY, CL, Cl,Cm, Cn
