#Inputs for the software


Flight Inputs:

    H_base = 0                            # The altitude in meters
    V_base = 70                           # The flight airspeed
    beta_base = 0 / 180 * math.pi         # The side slip angle
    gamma = 0                             # Flight path angle                                                                                              # previous condition  np.arctan(3 / 100)  # (3/100)/180*np.pi##math.atan(0/87.4)#/180*math.pi # 3% slope gradient # Best climb rate: 6.88m/s vertical @ 87.5m/s = 4.5°gamma, see http://www.atraircraft.com/products_app/media/pdf/Fiche_72-600_Juin-2014.pdf
    R = 000                               # in meters the turn radius   


Flight condition inputs


    Neng = 12                             # Number of engines
    inop_eng = 0                          # Number of inoperative engines
    FlapDefl = 0 * np.pi / 180            # in degree standard flap deflection. 
    g.VelFlap = 0                         # in m/s the maximum velocity at which flap are deployed
    alphastall = 11.7                     # that's for patterson




Geometry Inputs, class geometry

    # --- Mass ---
    x_cg = 11.75                            # Center of gravity position in meters from the tip of the aircraft
    m = 21500                               # Mass in Kg.


    # --- Geometry --- 
    
    WING
    S = 61.0                                # Wing surface m^2
    b = 27.05                               # Wingspan in m
    c = 2.324481                            # m Mean aerodynamic chord from OpenVSP, used as reference
    wingsweep = 0                           # Wing sweep in radians
    dihedral = 0                            # Dihedral angle in radians
    alpha_i = 4 / 180 * np.pi               # wing tilt angle, angle between reference line of fuselage and reference line of profile
    alpha_0 = -1.8/180*np.pi         # airfoil zero lift angle: from zero lift line to reference line. Negative means that airfoil lifts with 0 local angle of attack measured to reference line

    HT
    Sh = 11.13                              # Horizontal tail surface
    bh = 7.21                               # in m the HT wingspan
    c_ht = 1.54                             # Average chord of the horizontal tail
    self.zh = zh                            # Position of the horizontal tail on the vertical tail to the fuselage centre line
    Hor_tail_coef_vol = (Sh*lv) / (S*c)     # Volume coefficient of Horizontal tail
    it = -0.5 * np.pi/180                   # Horizontal tail tilt angle, with respect to the fuselage


    VT parameters
    self.Sv = Sv                            # Vertical tail surface
    self.bv = bv                            # Vertical tail wingspan
    self.Av = bv**2/Sv                      # Vertical tail aspect ratio
    self.bvl = bv+r                         # Vertical tailplane span extended to the fuselage center line
    fswept = 35/180*math.pi                 # sweep angle of Vertical Tail
    ftaper = 0.55                           # taper ratio of VT
    fAR = 1.57                              # Aspect ratio of VT (Wingspan/chord)


    FUSELAGE
    FusWidth = 2.82                         # Width of the fuselage
    self.r = r                              # Fuselage thickness at the section where the aerodynamic centre of vertical tail is
    self.zw = zw                            # wing position in fuselage. Height of wing root with respect to center of fuselage
    self.rf = rf                            # Fuselage max radius
    dfus=0,                                 # dfus=separation with the fuselage of first propeller NOT SURE IF RELATIVE
    dprop=0.1                               # Spacing between propellers       NOT SURE IF RELATIVE
    TipClearance=True  



    Propeller

    self.Dp                                 # Propeller diameter. Either is manually introduced or calculated with wingspan, N_eng, dprop, dfus and TipClearance
    self.Sp = self.Dp**2/4*math.pi          # Frontal surface propeller = pi*radio^2
    self.xp = self.Dp/2                     # Is the distance between propeller and leading edge
    self.step_y=self.Dp+dprop*self.Dp

    # --- Distances ---
    z_m = -0.443                            # vertical distance from center of gravity to propellers. Propellers are over Computed with OpenVSP
    z_h_w = 3.371                           # vertical distance from the horizontal tail to the propeller axis. Computed with OpenVSP
    lh = 14.2183                            # Horizontal distance between the aerodynamic centers of horizontal tail and wing (0.25 of their chord in root is enough) Computed with OpenVSP.
    lh2 = 12.09                             # Horizontal distance from the wing trailing edge to the horizontal tail quarter chord point. Computed with OpenVSP

    lv = 25.84 - x_cg                       # m Checked OpenVSP, distance from center of gravity to center of pressure of horizontal tail
    zf = 2                                  # z position of the MAC of the fin, in reality a suitable height
    lemac = 11.24                           # Distance from the tip to the leading edge of the MAC (here MAC matches with root chord)

    taudr = 0.30                            # Ruder efficiency factor see nicolosi paper and dela-vecchia thesis
    Var_xac_fus = -0.69                     # Variation in the aerodynamic centre (NP) between  congiguration on wing and wing + fuselage (without hor tail always) OpenVSP














    # Flap and aileron definition
    isflap = True                           # Condition for flap
    FlPosi = 0.05                           # Distance from start position of flap to center plane of aircraft divided by wingspan [0,0.5] 
    FlRatio = 0.75-FlPosi*2                 # Total flap length divided by wingspan
    FlChord = 0.3                           # Flap chord divided by total chord. (Not extandable flaps)

    isail = True                            # Condition for flap
    AilDiff = 0.5                           # Distance from start position of flap to center plane of aircraft divided by wingspan [0,0.5]
    AilPosi = 0.375                         # Distance from start position of flap to center plane of aircraft divided by wingspan [0,0.5]
    AilRatio = 0.125                        # Total flap length divided by wingspan
    AilChord = 0.3                          # Flap chord divided by total chord. (Not extandable flaps)               




    # Inertia terms are obtained from VSPaero for an homogeneous weight distribution
    Ix = 289873  # Kg/m^2
    Iy = 298442  # Kg/m^2
    Iz = 573579  # Kg/m^2
    Ixz = 0      # Kg/m^2


    # --- Power ---
    P_a = 2.051*10**6  # per engine
    P_var = P_a
    hp = 0  # rotor term
    prop_eff = 0.8 # Propeller efficiency
    ip = -1.6/180*np.pi  # propeller incidence angle with respect to zero lift line of the profile. Negative means propeller line is below zero lift line
    Pkeyword = 'Default'  # designate the propulsion model used to compute thrust






Coefficients


    # ---Unique coeff ---
    aht = 0.6131             # Horizontal effective tail lift coefficient. Effective means the influence of the rest of the aircraft is considered at 70m/s, alpha=0 (donwwash and tail dynamic pressure). Dimensioned with S.
    aht2 = 0.78082           # Horizontal tail lift coefficient, for the tail analysed alone. Dimensioned with S.
    Cm_alpha_wb = 1.173310   # Cm_alpha_wb from OpenVSP Aircraft without hor. tail
    cm_0_s = -0.0494         # = (0.2941)*Var_xac_fus/c  #zero lift pitching moment of the wing section at the propeller axis location. From the xlfr5 file, alpha = 0°
    # Cm_de = -8             # per rad, is constant for DECOL    You can use the one from STAB file, or this one



    # alpha = 0 coefficients from OpenVSP
    CD0T = 0.03383                   # from OpenVSP, parasitic zero lift drag. Total friction drag with both VT and HT     
    CD0T_wo_VT = 0.03112             # Friction drag without VT
    CL0 = 0.516688                   # Total CL0 including horizontal tail
    CL0_HT = -0.0284                 # Interpolated effective zero lift of horizontal tail (70 m/s). Effective means the influence of the rest of the aircraft is considered (donwwash and tail dynamic pressure)
    Cm0 = 0.035015                   #
    Cm0_wo_HT = -0.129536            # Cm0 of aircraft less horizontal tail


    # Drag polar without patterson. Interpolated from VSP v26, updated VSPAERO
    Cda_fl_0 = 1.1458               # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for no flaps
    Cdb_fl_0 = 0.1891
    Cdc_fl_0 = 0.026

    # with flaps down 15°, coefficients and drag polar
    Cd0_fl_15 = 0.030989             # extra drag when deflecting flaps 15°
    CL0_fl_15 = 0.500229             # extra lift when deflecting flaps 15°
    Cm0_fl_15 = -0.027008            # extra moment when deflecting flaps 15°

    Cda_fl_15 = 1.034                # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 15 °
    Cdb_fl_15 = 0.3506
    Cdc_fl_15 = 0.057

    # with flaps down 30°
    Cd0_fl_30 = 0.069758             # Extra lift when deflecting flaps 30°
    CL0_fl_30 = 0.862223             # Extra drag when deflecting flaps 30°
    Cm0_fl_30 = -0.046686            # Extra moment when deflecting flaps 30°

    Cda_fl_30 = 0.9197               # Coefficients for calculus of CD (CD=Cda * alpha ** 2 + Cdb * alpha + Cdc) for flaps = 30 °
    Cdb_fl_30 = 0.4424
    Cdc_fl_30 = 0.0957



    # Down-Wash parameters, measured from OpenVSP
    # No flaps
    eps0_flaps0 = 1.5 * np.pi/180         # Downwash at 0 angle of attack, no flaps
    deps_dalpha_flaps0 = 0.247            # Derivative of downwash with respect to alpha, no flaps 

    # 15° flaps
    eps0_flaps15 = 2.411                  # Downwash at 0 angle of attack flaps = 15°
    deps_dalpha_flaps15 = 0.2387          # Derivative of downwash with respect to alpha, for flaps = 15 °

    # 30° flaps
    eps0_flaps30 = 3.05 * np.pi/180       # Downwash at 0 angle of attack flaps = 30°
    deps_dalpha_flaps30 = 0.262           # Derivative of downwash with respect to alpha, for flaps = 30 °

    var_eps = 1.8                         # Parameter for inflow in slisptream. See Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof
    K_e = 1.35                            # Down wash factor, see Modeling the Propeller Slipstream Effect on Lift and Pitching Moment, Bouquet, Thijs; Vos, Roelof










Especial inputs

    # --- Propeller-Wing activation ---
    IsPropWing = True
    IsPropWingDrag = True

    
    # Airfoil characteristics
    Cd0_laminar = 0.009
    Cd0_turbulent = 0.009

    # Input file name
    Files = ['cldistribution', 'polar', 'flappolar', 'aileronpolar']  # best to replace the value
    alphaVSP = 5/180*np.pi

    PolarFlDeflDeg = 10              # Flap deflection for the naca3318fl+10 file used. File read in PattersonAugmented
    PolarAilDeflDeg = 10             # Aileron deflection for the naca3318fl+10 file used. File read in PattersonAugmented

    #Angles
    
    
    alpha_max = 15/180*np.pi         # Maximum angle for sections, after it, Jameson is used. Measured to fuselage reference. 
    alpha_max_fl = 10/180*np.pi      # Maximum angle for sections with flap, after it Jameson is used. Measured to fuselage reference.

    