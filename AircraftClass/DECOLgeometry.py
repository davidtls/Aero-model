
"""
Created on Thu Mar  7 14:58:23 2019

Geometry class for the DECOL

A few constants are defined here:
    N_eng
    inop_eng
    Vertical Tail correction terms

The rest are terms relative to DECOL geometry mass etc...

@author: e.nguyen-van
"""
import math
import sys
import numpy as np
from scipy.special import erf
from SeligPropellerReadPositiveT import Propu
from SeligPropellerRead import Propu
import os

class data:
    # all data from DECOL go here
    hangar = {}
    # shared data between class go here:


    # --- Mass ---
    x_cg = 0.713  # (m)
    m = 8.25  # Kg


    # --- Geometry --- 
    S = 0.5  # m^2
    b = 2  # m
    c = 0.25  # m    mean aerodynamic chord, used as reference
    lv = 1.809-x_cg  # m distance from center of gravity to center of pressure of horizontal tail
    zf = 0.165+0.085  # z position of the MAC of the fin, in reality a suitable height
    lemac = 0.610124  # Distance from the tip to the leading edge of the MAC (here MAC matches with root chord)
    fswept = 0/180*math.pi    # swept angle of vertical tail
    ftaper = 1.0                # Fin taper ratio, taper ratio of VT
    fAR = 1.8                   # Fin aspect ratio, aspect ratio of VT
    FusWidth = 0.170          # width of the fuselage
    bh = 0.611                # Horizontal tail wingspan
    Sh = 0.0877  # Horizontal tail surface
    Hor_tail_coef_vol = (Sh*lv) / (S*c)     # 0.7552  Volume coef of Horizontal tail
    it = 2.1/180*np.pi             # Horizontal tail tilt angle
    taudr = 0.24  # ruder efficiency factor see nicolosi paper and dela-vecchia thesis
    Var_xac_fus = -0.019914  # Variation in the aerodynamic centre. Compute by difference of the neutral points on OpenVSP between wing and wing + fuselage (without hor tail)

    wingsweep = 0  # radians,
    dihedral = 0


    # flap and aileron definition
    isflap = True
    FlPosi = 0.0425  # with respect to wingspan, the start position of flap [0,0.5]
    FlRatio = 0.2325*2  # the total flap length to wingspan ratio
    FlChord = 0.25  # with respect to local chord, flap chord
    isail = True
    AilDiff = 0.5  # ratio of aileron differential: deflection down is only half as important
    AilPosi = 0.275  # [0,0.5]
    AilRatio = 0.225  # One aileron should be consistent with b/2 >= AilRatio + AilPosi
    AilChord = 0.25




    # Inertia measured with Fernando
    Ix = 1.1  # Kg/m^2
    Iy = 1.2  # Kg/m^2
    Iz = 2.0  # Kg/m^2
    Ixz = 0   # Kg/m^2


    #--- Power ---
    P_a = 576      # in W equivalent of a twin engine. Limited to 10A per controller
    hp = 0  # rotor term
    prop_eff = 0.7       # Propulsion efficiency
    ip = (-2.41)/180*np.pi  # propeller incidence angle with respect to zero lift line
    Pkeyword = 'DefaultPatterson'  # designate the propulsion model used to compute thrust



    # --- Propeller-Wing activation ---
    IsPropWing=True
    IsPropWingDrag=True


    # --- Distances ---
    z_h_w = 0.025  # vertical distance from the horizontal tail 1/4 point to the propeller axis. Computed with OpenVSP. Positive if tail is over
    lh = 1.16575    # Horizontal distance between the aerodynamic centers of horizontal tail and wing (0.25 of their chord in root is enough) Computed with OpenVSP.
    lh2 = 0.942  # Horizontal distance from the wing trailing edge to the horizontal tail leading edge. Computed with OpenVSP
    K_e = 1.44   # Down wash factor, see Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof
    c_ht = 0.145  # Average chord of the horizontal tail
    var_eps = 1.5  # parameter for inflow in slisptream. See Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof
    cm_0_s = -0.0512  #    +  (0.2792)*Var_xac_fus/c     zero lift pitching moment of the wing section at the propeller axis location. From the xlfr5 file, alpha = 0°

    # ---Unique coeff ---
    aht = 0.5448
    aht2 = 0.7049       # Horizontal tail lift coefficient, for the tail analysed alone. Dimensioned with S.
    Cm_alpha_wb = 1.228506  # Cm_alpha_wb from OpenVSP Aircraft without hor. tail
    #  Cm_de = -2 # per rad, is constant for DECOL             You can use the one from STAB file, or this one
    Cm_alpha_fus = 0.015*180/np.pi


    # alpha=0 coeff

    # without flaps
    CD0T = 0.0636         # global one  extracted from flight not stab the file
    CDO_wo_VT = 0.0627
    CL0 = 0.44768
    CL0_HT = -0.002489  # Interpolated effective zero lift of horizontal tail (23.5 m/s). Effective means the influence of the rest of the aircraft is considered (donwwash and tail dynamic pressure)
    Cm0 = 0.004684
    Cm0_wo_HT = -0.056585  # Cm0 of aircraft less horizontal tail


    # Drag polar without patterson. Interpolated from VSP v26, updated VSPAERO
    Cda_fl_0 = 1.6372
    Cdb_fl_0 = 0.2934
    Cdc_fl_0 = 0.0637

    # with flaps down 15°
    Cd0_fl_15 = 0.026699          # Extra drag due to flap for prop wing interaction if no patterson used
    CL0_fl_15 = 0.27142           # Extra lift due to flaps when no patterson used
    Cm0_fl_15 = 0.096931          # Extra pitch moment due to flaps

    Cda_fl_15 = 1.4937      # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 15 °
    Cdb_fl_15 = 0.4327
    Cdc_fl_15 = 0.0905











    # Down-Wash parameters

    # No flaps
    eps0_flaps0 = 2.274 * np.pi/180   # Downwash at 0 angle of attack in no flaps configuration
    deps_dalpha_flaps0 = 0.281        # Derivative of downwash with respect to alpha, no flaps conf

    # 15° flaps
    eps0_flaps15 = 4.07 * np.pi/180   # Downwash at 0 angle of attack in 15° flaps configuration
    deps_dalpha_flaps15 = 0.281       # Derivative of downwash with respect to alpha, 15° flaps conf

 




    #Airfoil characteristics
    Cd0_laminar = 0.01252 
    Cd0_turbulent = 0.01252


    #wing tilt angle, angle between reference line of fuselaje and reference line of profile
    alpha_i = 3.2/180*np.pi
    
    #airfoil zero lift angle       angle between reference line of profile and zero lift line of airfoil. Negative means that airfoil lifts with 0 local angle of attack
    alpha_0 = -2.41/180*np.pi  # update with vsp file
    
    

    
    # Input file name
    Files = ['cldistribution', 'polar', 'flappolar', 'aileronpolar']  # best to replace the value
    alphaVSP = 0/180*np.pi
    PolarFlDeflDeg = 5
    PolarAilDeflDeg = 5










    #unique data go here:
    def CalcKf(self, bv, r):
        return 1.4685*(bv/(2*r))**(-0.143)
    
    def CalcKw(self, zw, rf):
        return -0.0131*(zw/rf)**2-0.0459*(zw/rf)+1.0026
        
    def CalcKh(self, zh, bvl, Sh, Sv):
        x = zh/bvl
        Khp = 0.906*x**2-0.8403*x+1.1475
        Khs = math.exp(0.0797*math.log(Sh/Sv)-0.0079)
        return 1+Khs*(Khp-1)
    
    def CalcKdr(self, Kf, Av):
#        Kdr=(1+((Kf-1)/2.2))*(1.33-0.09*Av) # for T-tail formula
        Kdr = 1.07*(1+(Kf-1)/2.2)  # for a body mounted H tail
        return Kdr
    
    def set_nofin(self, boolean):
        if type(boolean) == bool:
               self.nofin = boolean  # flag to use or not rudder
        else:
            sys.exit("Error type of 'nofin' isn't a boolean. Stopping")
            
    def loadAtmo(self):
#        filename='/home/e.nguyen-van/Documents/codesign-small-tail/Python/DECOL_Polar/si2py.txt'#'/home/e.nguyen-van/Documents/codesign-small-tail/Python/DECOL_Polar'
        sep = '\t'
        file = open("si2py.txt", 'r')
        vecname = file.readline()
        index = 0
        VariableList = []
        condition = True
        while condition:
            VariableList.append(vecname[index:vecname.index(sep,index)])
            if VariableList[-1] == 'kvisc':
                condition = False
            index=vecname.index(sep,index)+1
            
        units = file.readline() # skip units
        data = []
        VariableDic = {}
        for j in range(len(VariableList)):
            exec("VariableDic['"+VariableList[j]+"'] = []")  # initialize my variables
            
        for line in file:
            mylist = []
            element = ""
            for k in range(len(line)):
                if line[k] != '\t' and line[k] != '\n':
                    element = element+line[k]
                else:
                    mylist.append(element)
                    element = ""
            data.append(mylist)
        file.close()
        
        for i in range(len(data)):
            for k in range(len(data[0])-1):
                exec("VariableDic['"+VariableList[k]+"'].append({})".format(float(data[i][k])))
            
        return VariableDic

    def SetEngineNumber(self, N_eng, inop_eng, TipClearance, dfus, dprop):
        # Used to position the engine
        # must be recalled to change engine number
        # adjusts prop dia and engine position
        
        self.N_eng = N_eng  # number of engines
        self.inop = inop_eng  # number of inoperative engines
        self.dprop = dprop  # spacing between propeller.
        
        if N_eng % 2 != 0:
            sys.exit("Number of engine is not even")
        
        if N_eng > 2:
            #Compute Propeller Diameter
            """ Propeller diameter
            Patterson theory is limited at wing tip. Should let one prop diameter as margin
            Validation case require prop until wing tip.
            Inclusion of margin coefficients (dtip, dfus) and seperation factor(dprop)
            """
            
            if TipClearance == True:

                self.Dp = (self.b/2-self.FusWidth/2)/(N_eng/2+dfus+(N_eng/2-1)*dprop)
            else:
                #Put one engine at wing tip
                self.Dp = (self.b/2-self.FusWidth/2)/(N_eng/2+dfus-0.5+(N_eng/2-1)*dprop)
            
            self.Sp = self.Dp**2/4*math.pi
            self.x_offset = self.Dp/2     #This is distance from propeller to leading edge
            self.step_y = self.Dp+dprop*self.Dp

            
            if TipClearance == True:
                self.yp = np.arange(self.FusWidth/2+self.Dp*(dfus+0.5), self.b/2, self.step_y)
            else:
                self.yp = np.arange(self.FusWidth/2+self.Dp*(dfus+0.5), self.b/2+self.Dp/2, self.step_y)
                
            self.yp = np.append(-self.yp, self.yp)

            self.xp = np.full(self.N_eng, self.x_offset + (self.x_cg - self.lemac))
            self.yp = np.sort(self.yp)
            self.zp = np.full(self.N_eng, 0.075)  # 0.075 approximate vertical distance between engine and CG. Propellers is above. Computed with OpenVSP
            
        else:
            #default twin engine position
            self.step_y = 0.2
            self.Dp = 0.336  # 14" propeller
            self.Sp = self.Dp**2/4*math.pi
            self.x_offset = self.Dp/2       #This is distance from propeller to leading edge

            self.xp = np.full(self.N_eng, self.x_offset + (self.x_cg - self.lemac))
            self.yp = np.array([-self.step_y, self.step_y])
            self.zp = np.full(self.N_eng, 0.075)  # 0.075 approximate vertical distance between engine and CG. Propellers is above. Computed with OpenVSP
            
        return

    
    def __init__(self, VTsize, N_eng, inop_eng, FlapDefl, bv=0.329, r=0.057/2, zw=0.073, rf=0.085, zh=0, bvl=0.348, Sh=0.0877, Sv=0.06, TipClearance = True, dfus=0, dprop=0.1):
        self.VTsize=VTsize
        
        self.SetEngineNumber(N_eng, inop_eng, TipClearance, dfus, dprop)

        # See Nicolosi 2017, Ciliberti 2017, and Ciliberti thesis 2012 for more info
        self.Sv = Sv  # Vertical tail surface
        self.SvBase = Sv  # Vertical tail surface
        self.bv = bv  # Vertical tail wingspan
        self.r = r    # Fuselage thickness (radius) at the section where the aerodynamic centre of vertical tail is
        self.Av = bv**2/Sv  # Vertical tail aspect ratio
        self.Sh = Sh  # Horizontal tail surface
        self.zw = zw  # wing position in fuselage. Height of wing root with respect to center of fuselage
        self.rf = rf  # Fuselage max radius
        self.zh = zh  # Position of the horizontal tail on the vertical tail to the fuselage centre line
        self.bvl = bv+r  # Vertical tailplane span extended to the fuselage center line
        
        #Update drag
        self.Cdvt = 0.01374*Sv/self.S
        
        #Nicolosi csts
        self.Kf = self.CalcKf(bv, r)
        self.Kw = self.CalcKw(zw, rf)
        self.Kh = self.CalcKh(zh, bvl, Sh, Sv)
        self.Kdr = self.CalcKdr(self.Kf, self.Av)
        self.taudrBase = 0.018/Sv*0.8  #Rudder surface/fin surface * efficiency of rudder, rudder surface proportional to vtsize
        #Atmosphere
        self.AtmoDic = self.loadAtmo()

        #Aero coefficients for DragQuad
        self.FlapDefl = FlapDefl
        self.Cdo_fl, self.CL0_fl, self.Cm0_fl, self.Cda, self.Cdb, self.Cdc = self.AeroCoefs(FlapDefl)
        
        #Load propeller
        self.InitPropeller(8, 6)


    def AeroCoefs(self, FlapDefl):
         if FlapDefl == 0:
             self.Cd0_fl = 0
             self.CL0_fl = 0
             self.Cm0_fl = 0
             self.Cda = self.Cda_fl_0
             self.Cdb = self.Cdb_fl_0
             self.Cdc = self.Cdc_fl_0

             self.eps0 = self.eps0_flaps0
             self.deps_dalpha = self.deps_dalpha_flaps0
         elif FlapDefl == 15 * np.pi / 180:
             self.Cd0_fl = self.Cd0_fl_15
             self.CL0_fl = self.CL0_fl_15
             self.Cm0_fl = self.Cm0_fl_15
             self.Cda = self.Cda_fl_15
             self.Cdb = self.Cdb_fl_15
             self.Cdc = self.Cdc_fl_15

             self.eps0 = self.eps0_flaps15
             self.deps_dalpha = self.deps_dalpha_flaps15
         else:
             print("Chose an allowable value for flaps deflection, options are: No flaps (0°) or 15°")
             exit()
         return self.Cd0_fl, self.CL0_fl, self.Cm0_fl, self.Cda, self.Cdb, self.Cdc


    def InitPropeller(self, Dia, Pitch):                                       #Dia: Diameter in inches of propeller
        self.Propeller = Propu(Dia, Pitch)

        #update information about propeller
        self.Dp = Dia*0.0254
        self.Sp = np.pi*self.Dp**2/4

    # functions go here:
    def printVTsize(self):
        print('Asked vertical tail size of current class:')
        print(self.VTsize)
    
    def NicolosiCoef(self, MCoef, Mach):
        # function to compute Cy, Cl and Cn derivatives using VeDSC methods
        # replaces the coefficients in the matrix Mcoef by those computed by VeDSC
        # to call only once, it computes for every velocities/Mach number. Then linear interpol
        # Cy_beta, Cy_p = 0, Cy_r = 0, Cl_beta = 0, Cl_r = 0, Cn_beta= 0, Cn_p = 0, Cn_n = 0
        
        MVeDSC = np.copy(MCoef)  # create a copy of the coefficients matrix
        if self.nofin == False:
            # add a column to account for rudder
            dimZero = len(MVeDSC[:, 0])
            MVeDSC = np.hstack((MVeDSC, np.zeros((dimZero, 1))))
        
        K = self.Kf*self.Kh*self.Kw
        for i in range(len(Mach)):
            Av = self.bv**2/self.Sv
#            print(Av)
            cla = 2*math.pi # local lift curve slope coefficient of fin (perpendicular to leading edge) per rad
            eta = cla/(2*math.pi)
            av = -(Av * cla * math.cos(self.fswept) * 1/math.sqrt(1-Mach[i]**2*(math.cos(self.fswept))**2))/(Av*math.sqrt(1+4*eta**2/(Av/math.cos(self.fswept))**2)+2*eta*math.cos(self.fswept))# mach number formula
#            print(av)
            VeDSC_Coef = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
            VeDSC_Coef[0, 0] = K*av*self.Sv/self.S
            VeDSC_Coef[0, 1] = K*av*self.Sv/self.S*self.zf/self.b*2
            VeDSC_Coef[0, 2] = -K*av*self.Sv/self.S*2*self.lv/self.b #apparently x position of fin doesn't change
            VeDSC_Coef[0, 3] = -self.Kdr*av*self.taudr*self.Sv/self.S
            VeDSC_Coef[1, 0] = K*av*self.Sv/self.S*2*self.zf/self.b
            VeDSC_Coef[1, 1] = 0
            VeDSC_Coef[1, 2] = -K*av*self.Sv/self.S*self.zf/self.b*self.lv/self.b*2.0
            VeDSC_Coef[1, 3] = -self.Kdr*av*self.taudr*self.Sv/self.S*2*self.zf/self.b
            VeDSC_Coef[2, 0] = -K*av*self.lv/self.b*self.Sv/self.S
            VeDSC_Coef[2, 1] = -K*av*self.lv/self.b*self.Sv/self.S*self.zf/self.b*2.0
            VeDSC_Coef[2, 2] = K*av*self.lv/self.b*self.Sv/self.S*self.lv/self.b*2.0
            VeDSC_Coef[2, 3] = self.Kdr*av*self.taudr*self.lv/self.b*self.Sv/self.S
            # Coefficients are computed now access the right matrix and replace them
            VarPosi = (1, 2, 4)
            EffPosi = (1, 3, 5)
            NumEff = 6  # number of force equations
#            print(VeDSC_Coef[2,2])
            for kk in range(len(EffPosi)):
                # Replace rudder coefficient
                if self.nofin == False:
                    MVeDSC[EffPosi[kk]+i*NumEff, -1] = VeDSC_Coef[kk, -1]
                # Now coefficients are from the finless DECOL. Add new coefficients to the matrix
                for jj in range(len(VarPosi)):
                    if VeDSC_Coef[kk, jj] != 0:
                        MVeDSC[EffPosi[kk]+i*NumEff, VarPosi[jj]] = MVeDSC[EffPosi[kk]+i*NumEff, VarPosi[jj]]+VeDSC_Coef[kk, jj]
#            print(VeDSC_Coef)
            
        return MVeDSC
        
    
    def AdjustVT(self, VTNewSize = 1):
        # Adjust VT geometry parameters based on VTsize
        # Change only Kf, Kw and Kh, does not modify internal goem param
        
        # local variables
        Sv = self.SvBase
        bv = self.bv
        r = self.r
        Av = self.Av
        
        if VTNewSize != 1:
            Sv = Sv*VTNewSize
            taudr = self.taudrBase*VTNewSize
        else:
            Sv = Sv*self.VTsize
            taudr = self.taudrBase
        
        self.bv = (Av*Sv)**0.5
        self.zh = r+0.77*bv
        self.bvl = bv+r
        self.taudr = taudr
        
        # Compute new coef if necessary
        self.Sv = Sv
        self.Kf = self.CalcKf(self.bv, r)
        self.Kh = self.CalcKh(self.zh, self.bvl, self.Sh, Sv)
        
        
        return np.array([self.Kf, self.Kh])
    
    def AdjustVTcstSpan(self, VTNewSize=1):
        # Adjust VT geometry parameters based on VTsize
        # Based on Sv and cste span (increase Av)
        # Change only Kf, Kw and Kh
        
        # local variables
        Sv = self.SvBase
        bv = self.bv
        
        if VTNewSize != 1:
            Sv = Sv*VTNewSize
        else:
            Sv = Sv*self.VTsize
        
        self.Av = bv**2/Sv
        
        # Compute new coef if necessary
        self.Sv = Sv
        self.Kh = self.CalcKh(self.zh, self.bvl, self.Sh, Sv)
        
        return np.array([self.Kh])
    
    def GetAtmo(self, h=0):
        """
        Using the atmosphere model loaded before, it outputs [a_sound, rho] at
        the desired h=altitude. It doesn't perform interpolation.
        """
        Condition = h/500 is int
        if Condition:
            Indice = h//500+1
            
        else:
            if (h/500) < (h//500+500/2):
                Indice = int(h//500+1)
            else:
                Indice = int(h//500+2)
        
        results = np.array([self.AtmoDic['a'][Indice], self.AtmoDic['dens'][Indice]])
        return results
    
    def DefaultProp(self, dx, V):
        '''
        Compute thrust based on throttle levels dx
        Uses original model
        '''
        Thr = np.sum(dx)*2*self.P_var/(float(self.N_eng)*V)*self.prop_eff
        
        return Thr
    
    def DefaultPropPatterson(self, dx, V):
        #Same as defaultprop but returns a vector instead
        Thr = dx*2*self.P_var/(float(self.N_eng)*V)*self.prop_eff
        
        return Thr
    

    def DefaultTorque(self, dx, V):
        '''
        Compute torque based on default algo
        '''
        M_y = 0 # initialization
        M_y = -np.dot(dx, self.yp)
        # for i in range(n_eng):
        #     M_y=M_y+x[start+i]*g.step_y*(n_eng-i)
        # for i in range(n_eng):
        #     M_y=M_y-x[-i-1]*g.step_y*(n_eng-i)
            
        M_y = M_y*2*self.P_var/(float(self.N_eng)*V)*self.prop_eff
        
        return M_y
    







    def Thrust(self, dx, V, theta_i=np.array([None])):
        if self.Pkeyword == 'Default':
            return self.DefaultProp(dx, V)
        
        if self.Pkeyword == 'BaseFunction':
            if theta_i.any() is None:
                return self.BasisFunThrust(self.Theta_i, dx, V)
            else:
                return self.BasisFunThrust(theta_i, dx, V)
        if self.Pkeyword == 'Selig':
            return self.ThrustSelig(dx,V)

        if self.Pkeyword == 'DefaultPatterson':
            return self.DefaultPropPatterson(dx, V)
        
        else:
            return None



    def Torque(self, dx, V, theta_i=np.array([None])):
        if self.Pkeyword == 'Default':
            return self.DefaultTorque(dx, V)

        if self.Pkeyword == 'BaseFunction':
            if theta_i.any() is None:
                return self.BasisFunTorque(self.Theta_i, dx, V)
            else:
                return self.BasisFunTorque(theta_i, dx, V)

        if self.Pkeyword == 'Selig':
                return self.TorqueSelig(dx, V)

        if self.Pkeyword == 'DefaultPatterson':
            return self.DefaultTorque(dx, V)

        else:
            return None
        
    def DragModel(self, alpha, alpha_patter):
        # Added drag to account for induced drag due to flap deflection
        self.Cd_fl=(alpha_patter-alpha)*self.Cdalpha * self.FlRatio
        
        return None
    
    def ThrustSelig(self, dx, V):
        #Compute thrust of each engine based on propeller file from selig
        
        #interpol J first
        localJ = self.Propeller.DetermineJ(dx, V, self.P_var*2/self.N_eng)
        
        #get thrust
        Ct = self.Propeller.getCt(J=localJ)
        self.n = self.Propeller.getn(localJ, V)
        
        #get cp for tracking
        self.cp = self.Propeller.getCp(J=localJ)
        
        return Ct*1.225*self.n**2*self.Propeller.D**4
    
    def TorqueSelig(self, dx, V):
        
        Thr = self.ThrustSelig(dx, V)
        
        return -np.dot(self.yp, Thr)
