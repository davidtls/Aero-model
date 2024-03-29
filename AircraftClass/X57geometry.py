"""
Created on Saturday, October 1st 2022

Geometry class for the x57 Maxwell NASA


1) Flight cases


Flying condition explored:

First:
         Normal cruise. Just the two cruise engines working.
         Not blowing into the wing, so I guess no interaction considered. And Patterson is not applied to tip engines.
         Here we will be checking just OpenVSP.
         150 knots (77.1677 m/s) at 8000 ft (2438.4 m)
         Cruise CL = 0.75 at angle of attack 0°
         The predicted cruise wing drag coefficient is of 0.02191
         The total cruise is around CD_cruise = 0.0537249   ellos dicen 0.05423

         Propellers speed at speed cruise 2250 rpm


Second:
         Take-off/ Landing procedure with 30° deflection of flaps. 12 high lift propellers activated. Cruise propellers deactivated.
         There is a required lift coefficient of 3.95 for the stall speed of 58 knots (29.8378)

         The unblown maximum lift coefficient of the high-lift wing (with the 30° flap setting) is 2.439, which means
         that activating the high lift propellers should give the remaining  1.5. Results in the powered analysis publication
         show HLP can give up to 1.7.

3) Run the analysis with or without Vertical tail? Is VeDSC okay for this type of aircraft?

4) Define the flap and the control surfaces, in OpenVSP and XFLR5




Specifications (from https://www.nasa.gov/sites/default/files/atoms/files/x-57-litho-print-v4.pdf)
(based on Modification IV configuration)
Aircraft Weight – Approximately 3,000 pounds
Maximum Operational Altitude – 14,000 feet
Cruise Speed – 172 mph (at 8,000 feet)
Stall Speed – 58 knots (67 mph)
Batteries
• Lithium Ion
• 860 pounds
• 69.1 kilowatt hours (47 kilowatt hours usable)
Cruise Motors and Propellers (2)
• 60 kilowatts, continuous
• 72 kilowatts peak (at takeoff)
• Air-cooled
• Out-runner, 14-inch diameter
• 5-foot diameter propeller
• 117 pounds each, combined weight
High-Lift Motors and Propellers (12)
• 10.5 kilowatts
• Air-cooled
• In-runner, 6-inch diameter
• 5-blade, folding propeller
• 1.9 foot diameter propeller
• 15 pounds each, combined weight

The specially designed X-57 airfoil is tailored for a cruise lift coefficient of 0.90 and incorporates a 25% chord flap.
The flap design uses a single-pivot displaced hinge with a 30° maximum deflection.


In designing the X-57’s HLPS, several constraints were applied to simplify the design [15]. In particular, the high-lift
propellers’ longitudinal degrees of freedom were constrained to be equal: all propellers had the same offset relative to
the wing in the x and z directions†, and all nacelles were inclined to be aligned with the freestream in cruise


sceptor_cdr_day_2_package

               Pg 9
               Leading edge sweep, deg :1.9


               Pg 17
               DATA FROM OVERFLOW
               Torque N m : 22.4
               Power: 10.7 KW
               Thrust: 216.2
               Efficiency: 0.605


               Power each engine: 10.5 KW (peak) at the take-off

               The 10.3 kW power condition corresponds to the expected power required to produce 52.4 lbf of thrust
               at 58 KEAS which represents the high-lift propeller condition at the minimum steady, level flight speed


               The 4.9 kW condition corresponds to the expected power that the high-lift propeller requires according
               to the airspeed mode schedule at the reference approach speed of 75 KEAS, Vref goal of the X-57




Evaluation of Off-Nominal Performance and Reliability of a
Distributed Electric Propulsion Aircraft during Early Design, Pg 22

               2 x cruise motor 106.14 kg
               12 x cruise motor 81.65 kg
               Battery 390.08 kg
               Empennage 27.3 kg
               Fuselage 235.87 kg
               Landing gear 61.15 kg
               Wing 152.88 kg
               2 x pilot 170 kg
               Misc 135.7 kg
               Total 1360.77 kg



All data the presented here has been extracted from:

(VSP3 MODEL, PUBLICLY AVAILABLE) (1)
http://hangar.openvsp.org/vspfiles/414


(ALL GENERAL INFORMATION ABOUT X57)
https://www.nasa.gov/aeroresearch/X-57/technical/index.html



INFORMATION ABOUT PROPELLERS CHARACTERISTICS
               X-57 “Maxwell” High-Lift Propeller Testing and Model Development:
               https://ntrs.nasa.gov/api/citations/20210016834/downloads/LSAWT_HLP_Test_Aviation2021_Final0628.pdf

               Exploring the Effects of Installation Geometry in High-Lift Propeller Systems

               [19] Brandt, J. B., Deters, R. W., Ananda, G. K., and Selig, M. S., “UIUC Propeller Data Site,” http://mselig.
               ae.illinois.edu/props/propDB.html, 2017. URL http://m-selig.ae.illinois.edu/props/propDB.html, accessed:
               2017-05-31.

               [20] Patterson, M. D., and Borer, N. K., “Approach Considerations in Aircraft with High-Lift Propeller Systems,” 17th AIAA Aviation
               Technology, Integration, and Operations Conference, 2017.



VALIDATION OF AERODYNAMIC FORCES AND MOMENTS WITH RANS, FOR X-57 MOD-III
               - Computational Analysis of the External Aerodynamics of the Unpowered X-57 Mod-III Aircraft. Seung. Y. Yoo and Jared C. Duensing
               - Computational Analysis of a Wing Designed for the X-57 Distributed Electric Propulsion Aircraft
               Computational Simulations of Electric Propulsion Aircraft: the X-57 Maxwell (2019)
               - sceptor_cdr_day_2_package


@author: david.planas

"""
import math
import sys
import numpy as np




class data:
    # all data from x57 go here.
    hangar = {}
    # shared data between class go here:


    # --- Mass ---
    x_cg = 3.3560  # [m] (behind the tip)
    z_cg = 0.345948  # [m] (over the tip)
    m = 1360.77  # [Kg]


    # --- Geometry ---
    S = 6.196  # [m^2] Wing surface
    b = 9.642  # [m] Wingspan
    c = 0.6492451582  # [m] Mean aerodynamic chord used as reference. It is not the root chord, but is the mean chord (simple narrowing wing). Center of gravity 19% - 34%
    lv = 7.9049172 - x_cg  # [m] Checked OpenVSP, distance from center of gravity to center of pressure of horizontal tail
    zf = 0.5628  # [m] z position of the MAC of the fin, or of the 25% chord, with respect to center of gravity, in reality a suitable height. Here is positive if tail is over the cg.
    lemac = 3.179  # [m] Distance from the tip to the leading edge of the MAC (here MAC does not match with root chord, but is the chord whose chord = MAC, as is a simple narrowing wing) (3.104 from  aicraft's tip to leading edge of chord root) (chord of wing in root = 0.756m)
    fswept = 49/180*math.pi  # sweep angle of VT
    ftaper = 0.29  # Taper ratio of VT
    fAR = 1.34812  # Aspect ratio of VT
    FusWidth = 1.2192  # [m] In the location of the wing. Anyway this is important for placing the engines with the algorithm, not really if you place them manually
    bh = 3.1513  # [m] HT wingspan
    Sh = 2.4527  # [m^2] Horizontal tail surface
    Hor_tail_coef_vol = (Sh*lv) / (S*c)  # Volume coefficient of Horizontal tail
    it = 0 * np.pi/180  # [rad] Horizontal tail tilt angle.
    taudr = 0.5  # Ruder efficiency factor see nicolosi paper and dela-vecchia thesis. "A comprehensive review of vertical tail design"
    Var_xac_fus = -0.2678456   # [m] Variation in the aerodynamic centre. Compute by difference of the neutral points on OpenVSP between wing and wing + fuselage (without hor tail).
    # Fuselage moves forward the neutral point, so (Xnp(wing) - Xnp(wing+fuselage)) < 0

    wingsweep = 1.887*np.pi/180  # [rad] Sweep angle of the wing
    dihedral = 0*np.pi/180  # [rad] Dihedral angle of the wing


    # flap and aileron definition
    isflap = True
    FlPosi = 0.1265  # with respect to wingspan, the start position of flap [0,0.5]
    FlRatio = 0.5923  # the total flap length to wingspan ratio
    FlChord = 0.309  # with respect to local chord, flap chord
    isail = True
    AilDiff = 0.5
    AilPosi = 0.7198  # [0,0.5]
    AilRatio = 0.239/2  # One aileron.
    AilChord = 0.274




    # Inertia terms are obtained from "Flight dynamics and control assessment for differential thrust aircraft
    # in engine inoperative conditions including aero‑propulsive effects", TU DELFT, Maurice F.M. Hoogreef.
    Ix = 1325  # [Kg/m^2]
    Iy = 2161  # [Kg/m^2]
    Iz = 3193  # [Kg/m^2]
    Ixz = 0  # [Kg/m^2]



    # --- Distances ---
    z_h_w = -0.4494  # [m] vertical distance from the horizontal tail to the propeller axis. Computed with OpenVSP. Positive if tail is over.
    lh = 4.61  # [m] Horizontal distance between the aerodynamic centers of horizontal tail and wing (0.25 of their chord in root is enough) Computed with OpenVSP.
    lh2 = 4.057  # [m] Horizontal distance from the wing trailing edge to the horizontal tail leading edge. Computed with OpenVSP
    c_ht = 0.7782   # [m] Average chord of the horizontal tail



    cm_0_s = -0.1971  # zero lift pitching moment of the wing section (airfoil) at the propeller axis location. From the xlfr5 file, alpha = 0°



    # ---Unique coeff ---
    aht = 1.4026      # [1/rad] Horizontal effective tail lift coefficient. Effective means the influence of the rest of the aircraft is considered at 70m/s, alpha=0 (donwwash and tail dynamic pressure). Dimensioned with S.
    aht2 = 1.5578     # [1/rad] Horizontal tail lift coefficient, for the tail analysed alone. Dimensioned with S.
    Cm_alpha_wb = 0.0134 *180/np.pi  # [1/rad] Cm_alpha_wb from OpenVSP Aircraft without hor. tail.


    # alpha=0 coeff 77.67 m/s

    # without flaps
    CD0T = 0.0537249  # from analysis. OPENVSP30 gives 0.056028, parasitic zero lift drag. For flaps he says Cd0 = 0.0782 pATTERSON THESIS PG160
    CD0T_wo_VT = 0.0506759  # OpenVSP gives 0.003049
    CL0 = 0.75  # OpenVSP24 gives 0.758195    # Total CL0 including horizontal tail
    CL0_HT = -0.068519    # Interpolated effective zero lift of horizontal tail (70 m/s). Effective means the influence of the rest of the aircraft is considered (donwwash and tail dynamic pressure)
    Cm0 = 0.022087 + 0.25
    Cm0_wo_HT = -0.24786 + 0.25  # Cm0 of aircraft less horizontal tail


    # Drag polar without patterson. Interpolated from VSP v26, updated VSPAERO
    Cda_fl_0 = 1.3393            # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 30 °
    Cdb_fl_0 = 0.2845            # alpha in radians!!
    Cdc_fl_0 = 0.0461



    # with flaps down 30°
    Cd0_fl_30 = 0  # extra lift  #DEFINIR
    CL0_fl_30 = 0  # extra drag  # DEFINIR
    Cm0_fl_30 = 0  # extra moment  # DEFINIIR

    Cda_fl_30 = 0  # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 30 °   # DEFINIR
    Cdb_fl_30 = 0  # alpha in radians!!      # DEFINIR
    Cdc_fl_30 = 0  # DEFINIR


    # Down-Wash parameters

    # No flaps
    eps0_flaps0 = 2.230251 * np.pi/180   # Downwash at 0 angle of attack in no flaps configuration
    deps_dalpha_flaps0 = 0.1530          # Derivative of downwash with respect to alpha, no flaps conf



    # 30° flaps
    eps0_flaps30 = 4 * np.pi/180     # [rad] Downwash at 0 angle of attack with 30° flaps configuration. Computed with conventional flap, not Fowler
    deps_dalpha_flaps30 = 0.1530         # [rad/rad = °/° ] Derivative of downwash with respect to alpha, 30° flaps configuration. Computed with conventional flap, not Fowler







    # Airfoil characteristics
    Cd0_laminar = 0.0053
    Cd0_turbulent = 0.012


    # wing tilt angle, angle between reference line of fuselage and reference line of profile
    alpha_i = 0 / 180 * np.pi  # [rad]


    # airfoil zero lift angle: from zero lift line to reference line. Negative means that airfoil lifts with 0 local angle of attack measured to reference line
    alpha_0 = -7.25/180*np.pi  # [rad]
    # This alpha_0 is just used for determining g.alpha_max and g.alpha_max_fl. They both determine the stall of the
    # profiles. This is fine.



    # Input file name
    Files = ['cldistribution', 'polar', 'flappolar', 'aileronpolar']  # best to replace the value
    alphaVSP = 0/180*np.pi  # [rad] Angle of attack used during the analyses in STAB and FEM files in OpenVSP
    PolarFlDeflDeg = 30   # [degree] Flap deflection for the naca3318fl+10 file used. File read in PattersonAugmented
    PolarAilDeflDeg = 10  # [degree] Aileron deflection for the naca3318fl+10 file used. File read in PattersonAugmented

    path = 'X-57_STAB/'
    filenameNoFin = [path + 'Mach1.stab', path + 'Mach2.stab', path + 'Mach3.stab', path + 'Mach4.stab', path + 'Mach5.stab']


    PropPath = "./X-57_FEM/"
    PropFilenames = {'fem': [PropPath+"Mach1",
                             PropPath+"Mach2",
                             PropPath+"Mach3",
                             PropPath+"Mach4",
                             PropPath+"Mach5"],
                     'AirfoilPolar': PropPath+"Airfoil.txt",
                     'FlapPolar': PropPath+"Airfoil-flap.txt",
                     'AileronPolar': PropPath+"Airfoil-Aileron-10degree.txt"}









    def __init__(self, VTsize, N_eng, inop_eng, FlapDefl, HLP):
        self.VTsize = VTsize
        self.HLP = HLP

        self.SetEngineNumber(N_eng, inop_eng,HLP)

        # See Nicolosi 2017, Ciliberti 2017, and Ciliberti thesis 2012 for more info
        self.Sv = 3.9017  # Vertical tail surface
        self.SvBase = self.Sv  # Vertical tail surface
        self.bv = 1.4979  # Vertical tail wingspan
        self.r = 0.2232    # Fuselage thickness at the section where the aerodynamic centre of vertical tail is
        self.Av = self.bv**2/self.Sv  # Vertical tail aspect ratio
        self.Sh = 2.4527  # Horizontal tail surface
        self.zw = 0.9310  # wing position in fuselage. Height of wing root with respect to center of fuselage
        self.rf = 1.3732/2  # Fuselage max radius
        self.zh = 0  # Position of the horizontal tail on the vertical tail to the fuselage centre line
        self.bvl = self.bv+self.r  # Vertical tailplane span extended to the fuselage center line

        #Nicolosi csts
        self.Kf = self.CalcKf(self.bv, self.r)
        self.Kw = self.CalcKw(self.zw, self.rf)
        self.Kh = self.CalcKh(self.zh, self.bvl, self.Sh, self.Sv)
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
        x = zh/bvl
        Khp = 0.906*x**2-0.8403*x+1.1475
        Khs = math.exp(0.0797*math.log(Sh/Sv)-0.0079)
        return 1+Khs*(Khp-1)

    def CalcKdr(self, Kf, Av):
        # "A COMPREHENSIVE REVIEW OF VERTICAL TAIL DESIGN" Ciliberti 2017, Page 10
        Kdr=(1+((Kf-1)/2.2))*(1.33-0.09*Av) # for T-tail formula
        return Kdr

    def set_nofin(self, boolean):
        if type(boolean) == bool:
            self.nofin = boolean  # flag to use or not rudder
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








    def SetEngineNumber(self, N_eng, inop_eng,HLP):
        # Used to position the engine
        # must be recalled to change engine number
        # adjusts prop dia and engine position

        self.N_eng = N_eng  # number of engines
        self.inop = inop_eng  # number of inoperative engines
        self.HLP = HLP



        if HLP == True:


            self.Dp =np.full(self.N_eng,22.67 * 0.0254)
            self.Sp = self.Dp**2/4*math.pi


            self.xp = np.array([9.3, 11.6, 9.3, 11.6, 9, 10.5, 10.5, 9, 11.6, 9.1, 11.6, 9.1])*0.02547
            self.yp = np.array([-148.38, -125.7, -103.02, -80.34, -57.66, -34.98, 34.98, 57.66, 80.34, 103.02, 125.7, 148.38])*0.0254
            self.zp = np.full(self.N_eng, -0.454052)  # vertical distance from center of gravity to propellers. Computed with OpenVSP

            #self.x_offset = 10*0.0254
            self.x_offset = np.array([9.3, 11.6, 9.3, 11.6, 9, 10.5, 10.5, 9, 11.6, 9.1, 11.6, 9.1])*0.0254

            self.ip = np.full(self.N_eng, -7.25/180 * np.pi)
            # propeller incidence angle with respect to zero lift line of the profile. Negative means propeller line is below zero lift line.
            # 0 is the angle of the propellers with respect to horizontal (0°angle of attack). There is, however, profile twisting and you would still have to account for
            #the angle between the horizontal and the zero lift line of the propeller

        else:

            self.Dp = np.full(N_eng, 60 * 0.0254)
            self.Sp = self.Dp**2/4*math.pi

            self.xp = np.full(self.N_eng, 14.132*0.0254)  # Distance from the propeller to the tip of
            self.yp = np.array([-189.74, 189.74])*0.0254
            self.zp = np.full(self.N_eng, -0.5851)   # vertical distance from center of gravity to propellers. Computed with OpenVSP

            self.x_offset = np.full(N_eng, 0.3283)

            self.ip = np.full(self.N_eng, -7.25/180 *np.pi)


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
            print("Chose an allowable value for flaps deflection, options are: No flaps (0°) or 30°")
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
            cla = 2*math.pi  # local lift curve slope coefficient of fin (perpendicular to leading edge) per rad
            eta = cla/(2*math.pi)
            av = -(Av * cla * math.cos(self.fswept) * 1/math.sqrt(1-Mach[i]**2*(math.cos(self.fswept))**2))/(Av*math.sqrt(1+4*eta**2/(Av/math.cos(self.fswept))**2)+2*eta*math.cos(self.fswept))# mach number formula
            #            print(av)
            VeDSC_Coef = np.array([[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0]])
            VeDSC_Coef[0, 0] = K*av*self.Sv/self.S                                          # Cy_beta_VT
            VeDSC_Coef[0, 1] = K*av*self.Sv/self.S*self.zf/self.b*2                         # Cy_p_VT
            VeDSC_Coef[0, 2] = -K*av*self.Sv/self.S*2*self.lv/self.b                        # Cy_r_VT
            VeDSC_Coef[0, 3] = -self.Kdr*av*self.taudr*self.Sv/self.S                       # Cy_delta-r
            VeDSC_Coef[1, 0] = K*av*self.Sv/self.S*2*self.zf/self.b                         # # Cl_beta_VT
            VeDSC_Coef[1, 1] = 0                                                            # Cl_p_VT
            VeDSC_Coef[1, 2] = -K*av*self.Sv/self.S*self.zf/self.b*self.lv/self.b*2.0       # Cl_r_VT
            VeDSC_Coef[1, 3] = -self.Kdr*av*self.taudr*self.Sv/self.S*2*self.zf/self.b      # Cn_delta-r
            VeDSC_Coef[2, 0] = -K*av*self.lv/self.b*self.Sv/self.S                          # Cn_beta_VT
            VeDSC_Coef[2, 1] = -K*av*self.lv/self.b*self.Sv/self.S*self.zf/self.b*2.0       # Cn_p_VT
            VeDSC_Coef[2, 2] = K*av*self.lv/self.b*self.Sv/self.S*self.lv/self.b*2.0        # Cn_r_VT
            VeDSC_Coef[2, 3] = self.Kdr*av*self.taudr*self.lv/self.b*self.Sv/self.S         # Cn_delta-r
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









    def HLP_thrust2(self,V, atmo, n):
        # returns a vector

        J = V / (n * self.Dp)
        Thr = (atmo[1] * n**2 * self.Dp**4) * (-0.1084*J**2 - 0.1336*J + 0.3934)

        Cq = -0.0146*(J)**3 - 0.003*(J)**2 + 0.0025*(J) + 0.0452

        # Interpolation of High Lift Propellers Ct - J from XROTOR:
        # X-57 “Maxwell” High-Lift Propeller Testing and Model Development,  Fig 14

        return Thr, Cq

    def Get_n(self, n, args):

        dx = args[0]
        V = args[1]
        rho = args[2][1]
        Pmax = 10500*12

        Cq = -0.0146*(V/(n*self.Dp))**3 - 0.003*(V /(n*self.Dp))**2 + 0.0025*(V /(n*self.Dp)) + 0.0452
        Cp1 = 2*np.pi*Cq
        Cp2 = ((Pmax/self.N_eng)*dx)/(rho*n**3*self.Dp**5)

        return (Cp1 - Cp2)




    def Cruise_thrust(self, dx, V,atmo):
        # returns a vector

        Thr = ((2250 * 2*np.pi)/60) * 255 * dx * 0.8826 / V

        # The maximum torque is 255 N m . At cruise the torque is 177 N m .
        # Torque x angular speed is the power, max power is 60.082 KW.
        # At cruise (77.1677 m/s) the estimated drag is 954 N, so efficiency must be 0.8826
        # This calculation of thrust is supposed to be used around the cruise point:
        # V = 77.1677 m/s ,  if drag is well stimated in OpenVSP (CD_cruise = 0.0537249) then
        # that means dx = 0.6941 at cruise. Outside that condition it will not work adequately as the
        # efficiency depends on the advance parameter J.

        return Thr






    def HLP_thrust(self, dx, V,atmo):
        # returns a vector


        # Previous model
        # J = (V)(nD) ;  V = freestream speed, n = revs/s, D=diameter
        # Here Ct = T / (rho n^2 D^4)
        # Thr = (rho n^2 D^4) * (f(J)) * dx

        # This model is valid for the validation since they have done the conditions saying  Power = 0.39
        # What i undertand from that is that propellers are to the 39% of their nominal power/thrust, therefore dx = 0.39
        # This model works well since for dx = 1 and V = 28.3 m/s and sea level it predicts around 233 N
        # From "X-57 “Maxwell” High-Lift Propeller Testing and Model Development" Page 11
        # For example, the 10.3 kW power
        # condition corresponds to the expected power required to produce 52.4 lbf (233.08) of thrust at 58 KEAS which represents the
        # high-lift propeller condition at the minimum steady, level flight speed goal of the X-57.

        J = V / ((4800/60) * self.Dp)
        Thr = (atmo[1] * ((4800/60))**2 * self.Dp**4) * (-0.1084*J**2 - 0.1336*J + 0.3934)*dx



        # Interpolation of High Lift Propellers Ct - J from XROTOR:
        # X-57 “Maxwell” High-Lift Propeller Testing and Model Development,  Fig 14

        return Thr









    def Thrust(self, dx, V, atmo):

        """
        if self.HLP == True:
            return self.HLP_thrust(dx, V, atmo)
        """



        if self.HLP == True:
            from scipy.optimize import fsolve
            n = fsolve(self.Get_n, np.full(self.N_eng, 20), [dx, V, atmo])  # returns n in rev/s


            Mach_max = 0.4961
            wmax = (((atmo[0]*Mach_max)**2 - (V*(1+0.8))**2)**0.5)/(self.Dp*0.5)  # in rad/s
            nmax = wmax/(2*np.pi) # in rev/s
            for i in range(len(n)):
                if n[i]>nmax[i]:
                    n[i] = nmax[i]


            Thr, Cq = self.HLP_thrust2(V, atmo, n)

            Q = Cq * atmo[1] * n**2 * self.Dp**5
            Cp = 2*np.pi*Cq
            P = Cp * atmo[1] * n**3 * self.Dp**5
            Total_power = np.sum(P)

            for i in range(len(dx)):
                if dx[i] == 0:
                    Thr[i] = 0

            return Thr


        if self.HLP == False:
            return self.Cruise_thrust(dx, V, atmo)

        else:
            print('WARNING, HLP undefined in main')
            return None






    def DragModel(self, alpha, alpha_patter):
        # Added drag to account for induced drag due to flap deflection
        self.Cd_fl = (alpha_patter-alpha)*self.Cdalpha * self.FlRatio

        return None