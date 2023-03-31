"""
Created on Wednesday, November 17 2022

Geometry class for the D08

@author: david.planas

Conversion linear factor from an A320: 0.116935483

I need the airfoil and the model in Openvsp3
"""
import math
import sys
import numpy as np




class data:
    # all data from D08 go here.
    hangar = {}
    # shared data between class go here:


    # --- Mass ---
    x_cg = 1.9566  # [m] (behind the tip)
    m = 142.711  # [Kg] 161 han metido maz baterias


    # --- Geometry ---
    S =   # [m^2] Wing surface
    b = 3.9875  # [m] Wingspan
    c =   # [m] Mean aerodynamic chord used as reference. It is not the root chord, but is the mean chord (simple narrowing wing). Center of gravity 19% - 34%
    lv =  - x_cg  # [m] Checked OpenVSP, distance from center of gravity to center of pressure of horizontal tail
    zf =   # [m] z position of the MAC of the fin, or of the 25% chord, with respect to center of gravity, in reality a suitable height. Here is positive if tail is over the cg.
    lemac =   # [m] Distance from the tip to the leading edge of the MAC (here MAC does not match with root chord, but is the chord whose chord = MAC, as is a simple narrowing wing) (3.104 from  aicraft's tip to leading edge of chord root) (chord of wing in root = 0.756m)
    fswept = *math.pi  # sweep angle of VT
    ftaper =   # Taper ratio of VT
    fAR =   # Aspect ratio of VT
    FusWidth = 0.45531  # [m] In the location of the wing. Anyway this is important for placing the engines with the algorithm, not really if you place them manually
    bh = 1.46233  # [m] HT wingspan
    Sh =   # [m^2] Horizontal tail surface
    Hor_tail_coef_vol = (Sh*lv) / (S*c)  # Volume coefficient of Horizontal tail
    it =  * np.pi/180  # [rad] Horizontal tail tilt angle.
    taudr =   # Ruder efficiency factor see nicolosi paper and dela-vecchia thesis. "A comprehensive review of vertical tail design"
    Var_xac_fus =    # [m] Variation in the aerodynamic centre. Compute by difference of the neutral points on OpenVSP between wing and wing + fuselage (without hor tail).
    # Fuselage moves forward the neutral point, so (Xnp(wing) - Xnp(wing+fuselage)) < 0

    wingsweep = 25*np.pi/180  # [rad] Sweep angle of the wing #27.45° I measure
    dihedral = 10.5*np.pi/180  # [rad] Dihedral angle of the wing


    # flap and aileron definition
    isflap = True
    FlPosi =   # with respect to wingspan, the start position of flap [0,0.5]
    FlRatio = -FlPosi*2  # the total flap length to wingspan ratio
    FlChord =   # with respect to local chord, flap chord
    isail = True
    AilDiff =
    AilPosi =   # [0,0.5]
    AilRatio =   # One aileron.
    AilChord =




    # Inertia terms are obtained from VSPaero for an homogeneous weight distribution
    Ix = 22.9245  # [Kg/m^2]
    Iy = 102.3728  # [Kg/m^2]
    Iz = 119.3115  # [Kg/m^2]
    Ixz = 6.16     # [Kg/m^2]


    # --- Power ---
    P_a =   # per engine
    P_var = P_a
    hp =   # rotor term
    prop_eff =
    ip = -3/180*np.pi  + alpha_zero # propeller incidence angle with respect to zero lift line of the profile. Negative means propeller line is below zero lift line
    Pkeyword = 'Default'  # designate the propulsion model used to compute thrust


    # --- Propeller-Wing activation ---
    IsPropWing = False
    IsPropWingDrag = False


    # --- Distances ---
    z_m =    # vertical distance from center of gravity to propellers. Propellers are over Computed with OpenVSP
    z_h_w =   # vertical distance from the horizontal tail to the propeller axis. Computed with OpenVSP
    lh =    # Horizontal distance between the aerodynamic centers of horizontal tail and wing (0.25 of their chord in root is enough) Computed with OpenVSP.
    lh2 =     # Horizontal distance from the wing trailing edge to the horizontal tail quarter chord point. Computed with OpenVSP
    c_ht =     # Average chord of the horizontal tail
    cm_0_s =  # +  (0.2941)*Var_xac_fus/c  #zero lift pitching moment of the wing section at the propeller axis location. From the xlfr5 file, alpha = 0°

    # ---Unique coeff ---
    aht =         # Horizontal effective tail lift coefficient. Effective means the influence of the rest of the aircraft is considered at 70m/s, alpha=0 (donwwash and tail dynamic pressure). Dimensioned with S.
    aht2 =        # Horizontal tail lift coefficient, for the tail analysed alone. Dimensioned with S.
    Cm_alpha_wb =  # Cm_alpha_wb from OpenVSP Aircraft without hor. tail



    # alpha=0 coeff

    # without flaps
    CD0T =   # from OpenVSP, parasitic zero lift drag      THESIS HAMBURG 0.027403
    CD0T_wo_VT =
    CL0 =       # Total CL0 including horizontal tail
    CL0_HT =    # Interpolated effective zero lift of horizontal tail (70 m/s). Effective means the influence of the rest of the aircraft is considered (donwwash and tail dynamic pressure)
    Cm0 =
    Cm0_wo_HT =     # Cm0 of aircraft less horizontal tail


    # Drag polar without patterson. Interpolated from VSP v26, updated VSPAERO
    Cda_fl_0 =
    Cdb_fl_0 =
    Cdc_fl_0 =


    # with flaps down 30°
    Cd0_fl_30 =     # extra lift
    CL0_fl_30 =     # extra drag
    Cm0_fl_30 =    # extra moment

    Cda_fl_30 =       # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 30 °
    Cdb_fl_30 =
    Cdc_fl_30 =


    # Down-Wash parameters

    # No flaps
    eps0_flaps0 =  * np.pi/180   # Downwash at 0 angle of attack in no flaps configuration
    deps_dalpha_flaps0 =       # Derivative of downwash with respect to alpha, no flaps conf

    # 30° flaps
    eps0_flaps30 = * np.pi/180
    deps_dalpha_flaps30 =


    # Airfoil characteristics
    Cd0_laminar = 0.009
    Cd0_turbulent = 0.009


    # wing tilt angle, angle between reference line of fuselage and reference line of profile
    alpha_i =  / 180 * np.pi


    # airfoil zero lift angle: from zero lift line to reference line. Negative means that airfoil lifts with 0 local angle of attack measured to reference line
    alpha_0 = /180*np.pi



    # Input file name
    Files = ['cldistribution', 'polar', 'flappolar', 'aileronpolar']  # best to replace the value
    alphaVSP = /180*np.pi
    PolarFlDeflDeg =    # Flap deflection for the naca3318fl+10 file used. File read in PattersonAugmented
    PolarAilDeflDeg =   # Aileron deflection for the naca3318fl+10 file used. File read in PattersonAugmented
    alpha_max = /180*np.pi
    alpha_max_fl = /180*np.pi


    path = 'D08/'
   filenameNoFin = [path + 'Mach1.stab', path + 'Mach2.stab', path + 'Mach3.stab', path + 'Mach4.stab', path + 'Mach5.stab']


    PropPath = "./D08/"
    PropFilenames = {'fem': [PropPath+"Mach1",
                         PropPath+"Mach2",
                         PropPath+"Mach3",
                         PropPath+"Mach4",
                         PropPath+"Mach5"],
                 'AirfoilPolar': PropPath+"Airfoil.txt",
                 'FlapPolar': PropPath+"Airfoil-flap.txt",
                 'AileronPolar': PropPath+"Airfoil-Aileron-10degree.txt"}











    def __init__(self, VTsize, N_eng, inop_eng, FlapDefl):
        self.VTsize = VTsize

        self.SetEngineNumber(N_eng, inop_eng)

        # See Nicolosi 2017, Ciliberti 2017, and Ciliberti thesis 2012 for more info
        self.Sv =   # Vertical tail surface
        self.SvBase = self.Sv  # Vertical tail surface
        self.bv =   # Vertical tail wingspan
        self.r =     # Fuselage thickness at the section where the aerodynamic centre of vertical tail is
        self.Av = self.bv**2/self.Sv  # Vertical tail aspect ratio
        self.Sh =   # Horizontal tail surface
        self.zw =   # wing position in fuselage. Height of wing root with respect to center of fuselage
        self.rf = 0.45531/2  # Fuselage max radius
        self.zh =   # Position of the horizontal tail on the vertical tail to the fuselage centre line
        self.bvl = self.bv+self.r  # Vertical tailplane span extended to the fuselage center line

        #Nicolosi csts
        self.Kf = self.CalcKf(bv, r)
        self.Kw = self.CalcKw(zw, rf)
        self.Kh = self.CalcKh(zh, bvl, Sh, Sv)
        self.Kdr = self.CalcKdr(self.Kf, self.Av)
        self.taudr = 0.30
        #Atmosphere
        self.AtmoDic = self.loadAtmo()

        #Aero coefficients for DragQuad
        self.FlapDefl = FlapDefl
        self.Cdo_fl, self.CL0_fl, self.Cm0_fl, self.Cda, self.Cdb, self.Cdc = self.AeroCoefs(FlapDefl)




        #unique data go here:
    def CalcKf(self, bv, r):
        # "A new vertical tailplane design procedure through cfd" Ciliberti thesis 2012, Page 111
        return 1.4685*(bv/(2*r))**(-0.143)  # Fuselage Correction factor.

    def CalcKw(self, zw, rf):
        # "A new vertical tailplane design procedure through cfd" Ciliberti thesis 2012, Page 112
        return -0.0131*(zw/rf)**2-0.0459*(zw/rf)+1.0026  # Wing correction factor

    def CalcKh(self, zh, bvl, Sh, Sv):
        # Horizontal tailplane correction factor
        # "A new vertical tailplane design procedure through cfd" Ciliberti thesis 2012, Page 113-114
        x=zh/bvl
        Khp=0.906*x**2-0.8403*x+1.1475
        Khs=math.exp(0.0797*math.log(Sh/Sv)-0.0079)
        return 1+Khs*(Khp-1)

    def CalcKdr(self, Kf, Av):
        # "A COMPREHENSIVE REVIEW OF VERTICAL TAIL DESIGN" Ciliberti 2017, Page 10
        Kdr=(1+((Kf-1)/2.2))*(1.33-0.09*Av) # for T-tail formula
        return Kdr

    def set_nofin(self, boolean):
        if type(boolean) == bool:
            self.nofin = boolean # flag to use or not rudder
        else:
            sys.exit("Error type of 'nofin' isn't a boolean. Stopping")

    def loadAtmo(self):
        filename = 'si2py.txt'
        sep = '\t'
        file = open(filename, 'r')
        vecname = file.readline()
        index = 0
        VariableList = []
        condition = True
        while condition:
            VariableList.append(vecname[index:vecname.index(sep, index)])
            if VariableList[-1] == 'kvisc':
                condition = False
            index = vecname.index(sep, index)+1

        units = file.readline()  # skip units
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







    def SetEngineNumber(self, N_eng, inop_eng):
        # Used to position the engine
        # must be recalled to change engine number
        # adjusts prop dia and engine position

        self.N_eng = N_eng  # number of engines
        self.inop = inop_eng  # number of inoperative engines

        self.PosiEng = np.array([-1.9938, -1.107, -0.66, 0.66, 1.107, 1.9938])
        self.Dp = np.array([0.3188, 0.4064, 0.4064, 0.4064, 0.4064, 0.3188])  # This must be a constant, not a vector! Otherwise code must be changed
        self.Sp = self.Dp**2/4*math.pi

        self.xp = np.array([ , ,  , , , ])
        self.yp = np.array([-1.9938, -1.107, -0.66, 0.66, 1.107, 1.9938])
        self.zp = np.array([ ,  ,  ,   ,   , ])  # vertical distance from center of gravity to propellers. Computed with OpenVSP

        self.x_offset = np.array([0.3107, 0.3166, 0.3166, 0.3166, 0.3166,0.3107 ])

        return




    # functions go here:
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
        elif FlapDefl == 30 * np.pi / 180:
            self.Cd0_fl = self.Cd0_fl_30
            self.CL0_fl = self.CL0_fl_30
            self.Cm0_fl = self.Cm0_fl_30
            self.Cda = self.Cda_fl_30
            self.Cdb = self.Cdb_fl_30
            self.Cdc = self.Cdc_fl_30

            self.eps0 = self.eps0_flaps30
            self.deps_dalpha = self.deps_dalpha_flaps30
        else:
            print("Chose an allowable value for flaps deflection, options are: No flaps, 15° or 30°")
            exit()
        return self.Cd0_fl, self.CL0_fl, self.Cm0_fl, self.Cda, self.Cdb, self.Cdc

    def printVTsize(self):
        print('Asked vertical tail size of current class:')
        print(self.VTsize)

    def NicolosiCoef(self, MCoef, Mach):
        # function to compute Cy, Cl and Cn derivatives using VeDSC methods
        # replaces the coefficients in the matrix Mcoef by those computed by VeDSC
        # to call only once, it computes for every velocities/Mach number. Then linear interpol
        # Cy_beta, Cy_p = 0, Cy_r = 0, Cl_beta = 0, Cl_r = 0, Cn_beta= 0, Cn_p = 0, Cn_n = 0

        MVeDSC=np.copy(MCoef) # create a copy of the coefficients matrix
        if self.nofin == False:
            # add a column to account for rudder
            dimZero = len(MVeDSC[:, 0])
            MVeDSC = np.hstack((MVeDSC, np.zeros((dimZero, 1))))

        K = self.Kf*self.Kh*self.Kw
        #        print('New VT efficiency = {0:0.4f}'.format(K))
        for i in range(len(Mach)):
            Av = self.bv**2/self.Sv
            #            print(Av)
            cla = 2*math.pi # local lift curve slope coefficient of fin (perpendicular to leading edge) per rad
            eta = cla/(2*math.pi)
            av = -(Av * cla * math.cos(self.fswept) * 1/math.sqrt(1-Mach[i]**2*(math.cos(self.fswept))**2))/(Av*math.sqrt(1+4*eta**2/(Av/math.cos(self.fswept))**2)+2*eta*math.cos(self.fswept))# mach number formula
            #            print(av)
            VeDSC_Coef = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0]])
            VeDSC_Coef[0, 0] = K*av*self.Sv/self.S
            VeDSC_Coef[0, 1] = K*av*self.Sv/self.S*self.zf/self.b*2
            VeDSC_Coef[0, 2] = -K*av*self.Sv/self.S*2*self.lv/self.b
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



            for kk in range(len(EffPosi)):
                #Manually change rudder coefficient by simple proportionality
                #MVeDSC[EffPosi[kk]+i*NumEff,-1]=MVeDSC[EffPosi[kk]+i*NumEff,-1]*self.VTsize
                # Replace rudder coefficient
                if self.nofin == False:
                    MVeDSC[EffPosi[kk]+i*NumEff, -1] = VeDSC_Coef[kk, -1]                  #Adds CY_delta_r   ,  Cl_delta_r  ,  Cn_delta_r  (there were 0s before)

                # Now coefficients are from the finless ATR. Add new coefficients to the matrix
                for jj in range(len(VarPosi)):
                    if VeDSC_Coef[kk, jj] != 0:
                        MVeDSC[EffPosi[kk]+i*NumEff, VarPosi[jj]] = MVeDSC[EffPosi[kk]+i*NumEff, VarPosi[jj]]+VeDSC_Coef[kk, jj]      #adds a term to CY_beta  CY_p CY_r
                        # Cl_beta  Cl_p  Cl_r   Cn_beta  Cn_p  Cn_r
        #            print(VeDSC_Coef)

        return MVeDSC


    def AdjustVT(self, VTNewSize=1):
        # Adjust VT geometry parameters based on VTsize
        # Change only Kf, Kw and Kh, does not modify internal goem param

        # local variables
        Sv = self.SvBase
        bv = self.bv
        r = self.r
        Av = self.Av

        if VTNewSize != 1:
            Sv = Sv*VTNewSize
        else:
            Sv = Sv*self.VTsize

        self.bv = (Av*Sv)**0.5
        self.zh = r+0.77*bv
        self.bvl = bv+r

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
        print('New VT aspect ratio Av={0:0.3f}'.format(self.Av))

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



    def ThrustCalculus(self, dx, V,atmo):
         # returns a vector

        # Some windmill is considered, for a very small range of J, after Thrust is made 0 artificially

        # dx = [engine left wing, ..., engine right wing]
        # V is a vector with size len(dx)

        rho =atmo[1]
        Thr = np.zeros(len(dx))

        for i in range(len(dx)):

             if i == 0 or i == 5:

                 # OUTER ENGINES
                 n = (5311.1*dx + 1442.8)/60 # [rps]
                 J = V[i]/(n*self.Dp)

                 if J <0.7125:
                     Ct = 0.35149
                     Thr[i] = Ct * (rho* (n)**2 * (self.Dp)**4)
                 if J>2.4:
                     Thr[i] = 0 # NO WINDMILLING CONSIDERED
                 else:
                     Ct = -0.1709*(J)**6 + 1.5416*(J)**5 - 5.4012*(J)**4 + 9.2325*(J)**3 - 8.1116*(J)**2 + 3.4855*(J) -0.2252
                     Thr[i] = Ct * (rho* (n)**2 * (self.Dp)**4)

             else:

                # INNER ENGINES
                n = (4723*dx + 638.48)/60 # [rps]
                J = V[i]/(n*self.Dp)

                if J <0.7125:
                    Ct = 0.35
                    Thr[i] = Ct * (rho* (n)**2 * (self.Dp)**4)

                if J>2.403:
                    Thr[i] = 0 # NO WINDMILLING CONSIDERED
                else:
                    Ct = 0.1476*(J)**6 - 1.4143*(J)**5 + 5.7064*(J)**4 - 12.352*(J)**3 + 14.717*(J)**2 - 8.9508*(J) + 2.4939
                    Thr[i] = Ct * (rho* (n)**2 * (self.Dp)**4)



        return Thr





    def Thrust(self, dx, V, ):

        if self.Pkeyword == 'Default':
            return self.ThrustCalculus(dx, V)
        else:
            print('WARNING, Pkeyword undefined in func "Thrust"')
            return None





    def DragModel(self, alpha, alpha_patter):
        # Added drag to account for induced drag due to flap deflection
        self.Cd_fl = (alpha_patter-alpha)*self.Cdalpha * self.FlRatio

        return None