from LaunchVehicle import LaunchVehicle
import math
from math import sqrt, pow, pi, cos, sin, tan, acos, asin, atan, radians, degrees, exp, log

class Step(LaunchVehicle):
    

    def __init__(self, LV, r, step_num, TW, engine, dome_shape, propellants, pressurant, insulation_dict, tank_material, num_engines, num_gimballed_engines, t_start, fairing_material, fairing_shape, num_boosters, pressure_fed):
        pass
        self.LV = LV
        self.r = r
        self.step_num = step_num
        self.TW = TW
        self.engine = engine
        self.dome_shape = dome_shape
        self.propellants = propellants
        self.pressurant = pressurant
        self.insulation_dict = insulation_dict
        self.tank_material = tank_material
        self.tank_ox_material = self.tank_material[0]
        self.tank_f_material = self.tank_material[1]
        self.num_engines = num_engines
        self.num_gimballed_engines = num_gimballed_engines
        self.t_start = t_start
        self.num_boosters = num_boosters
        self.pressure_fed = pressure_fed
        if ((self.LV.name == "Zephyr-1") | (self.LV.name == "Zephyr-2")):
            self.p_press = 500 * 0.101325 # 500 atmospheres in MPa in pressurant tank
        else:
            self.p_press = 20 * 0.101325 # 20 atmospheres in MPa in pressure tank
        if self.num_boosters == 0:
            self.parallel = False
        else:
            self.parallel = True
        
        if step_num == LV.num_steps:
            self.fairing_material = fairing_material
            self.fairing_shape = fairing_shape
        else:
            self.m_PL = 0
            self.fairing_material = 0
            self.fairing_shape = 0
        if self.parallel:
            self.fairing_material = fairing_material
            self.fairing_shape = fairing_shape

    def print(self):
        print("This is step " + str(self.step_num) + " of " + self.LV.name)
        print("The radius is " + str(self.r))
        print("The Thrust-to-Weight is "+ str(self.TW))
        print("The engine is based on the " + self.engine)
        print("The dome shape is " + self.dome_shape)
        print("The propellant combination is " + self.propellants)
        print("The body material is " + self.LV.body_material)
        print("The tank material is " + self.tank_material)
        print("The number of engines is " + str(self.num_engines))
        if self.step_num == self.LV.num_steps:
            print("The fairing material is " + self.fairing_material)
            print("The fairing shape is " + "'" + self.fairing_shape + "'")

    def sizeStep(self):
    ## Assumptions
    # Elliptical-2 dome has AR = 2 and Elliptical-sqr2 has AR = sqrt(2)
    # thickness of fairing is 0.005 m
    # Bulkheads are uncommon
    # volume of single cubesat is 0.001 m^3
    # 12 cubesats for mission 1 and 21 for mission 2
    # payload volume margin is 1.5: accounts for air gaps between cubeSats and
    # fairing wall
    #STEPSIZING Summary of this function goes here
    # The inputs of this function are:
        # the step numbers, countng from the 1st step
        # the gross lift-mass in kg => m_gross in kg
        # TW - Thrust to weight ratio for the stage
        # m_prop =  propellant mass(kg)
        # the engine used to size engine_params => engine = "Raptor" or
        # "Merlin"
        # type of propellant tank dome => dome = "Elliptical-2" (hemispheres) or "Elliptical-sqr2" or
        # "Circular" (hemispheres) or "Spherical"
        # type of propellants => propellants = "Methalox" or "Kerolox"
        # material of the tanks for shrinkage => right now the only option
        # is => tank_material = "Aluminum 6061"
        # number of engines for this step => num_engines = 1, 2, 3, 4, etc...
        # the material of the fairing is fairing_material = "Aluminum 6061",
        # for example
        # the shape of the fairing is fairing_shape = "Cone", "Ogive", etc.
        # the payload mass is m_PL, INCLUDING the PL attach fitting
        #self.step_params = {'Radius': None, 'Exit diameter': None, 'Nozzle length': None, 'Thrust structure': None, 'Aft skirt': None, 'Dome_f': None, 'Cyl_f': None, 'Intertank': None, 'Dome_ox': None, 'Cyl_ox': None, 'Fairing': None, 'Fwd skirt': None}
        # [length, SA, Vol]
        self.aft_skirt = []
        self.dome_f = []
        self.cyl_f = []
        self.intertank = []
        self.dome_ox = []
        self.cyl_ox = []
        self.srm_casing = []
        self.press_tank = []
        #if self.step_num == self.LV.num_steps:
        #    self.fairing = []
        #    self.fwd_skirt = []
        #if self.parallel:
        #    self.fairing = []
        self.fairing = []
        self.fwd_skirt = []

        # SET PARAMETERS
        pi = 3.1415926535897932#3846264338327950288
        g = 9.80665  # gravitational acceleration (m/s^2)
        p_inf = 0.101325  # atmospheric pressure (MPa)
        # if ((self.LV.name == "Zephyr-1") | (self.LV.name == "Zephyr-2")):
        #     p_tank = 30*p_inf # Pressure-fed raptor engines require tank pressure ~= chamber pressure
        # else:
        p_tank = 3 * p_inf # propellant tank pressure (MPa)
        rm_temp = 293  # room temperature (K)
        R_uni = 8314.4 # universal gas constant (J/kmol K)
        if self.engine == "Raptor":
            Isp = 330  # specific impulse (s) ACTUAL
            #Isp = 365 # specific impulse (s) ON EXCEL 
            Isp_vac = 380  # vacuum specific impulse (s)
            self.p_c = 30*pow(10,6)  # chamber pressure (Pa) ACTUAL
            self.epsilon = 45  # expansion ratio
        elif self.engine == "Merlin":
            Isp = 282  # specific impulse (s)
            Isp_vac = 311  # vacuum specific impulse (s)
            self.p_c = 9.7*pow(10,6)  # chamber pressure (Pa)
            self.epsilon = 22  # expansion ratio
        elif self.engine == "SRBNozzle": # Content is UNEDITED. NEEDS EDIT
            Isp = 265
            Isp_vac = 285
            self.p_c = 9.7*pow(10,6)
            self.epsilon = 22

        if self.dome_shape == "Elliptical-2":
            AR = 2  # aspect ratio ( > 1 means elliptical hemisphere, = 1 means spherical hemisphere) (1)
        elif self.dome_shape == "Elliptical-sqr2":
            AR = sqrt(2) 
        elif self.dome_shape == "Circular":
            AR = 1 
        e = sqrt(1 - pow(1/AR,2))  # eccentricity (1)

        self.setInsulation()

        if self.propellants == "Kerolox":
            fuel_temp = 293 
            ox_temp = 80 
            self.propulsion = 'Liquid'
            self.rho_f = 820  # fuel densitity (kg/m^3)
            self.rho_ox = 1140   # oxidizer density (kg/m^3)
            OF = 2.34  # O/F mixture ratio (1)
            FO = 1/OF  # F/O mixture ratio (1)
            self.fuel_frac = FO/(1 + FO)  # fraction of fuel mass (1)
        elif self.propellants == "Methalox":
            fuel_temp = 100  # current temperature of fuel(K)
            ox_temp = 80  # boiling temp of ox (K)
            self.propulsion = 'Liquid'
            self.rho_f = 424  # fuel densitity (kg/m^3)
            self.rho_ox = 1140   # oxidizer density (kg/m^3)
            OF = 3.55  # O/F mixture ratio (1)
            FO = 1/OF  # F/O mixture ratio (1)
            self.fuel_frac = FO/(1 + FO)  # fraction of fuel mass (1)
        elif self.propellants == "AP-Al-HTPB":
            prop_temp = rm_temp
            self.propulsion = 'Solid'
            self.rho_prop = 1854.42  # propellant density (kg/m^3) https://www.researchgate.net/figure/Formulation-of-propellant-A-The-density-of-cured-HTPB-is-920-kg-m-3_tbl1_279252447

        if self.pressurant == "Helium":
            Mol_mass = 4.0026 # molecular mass (kg/kmol)
        elif self.pressurant == "Nitrogen":
            Mol_mass = 14.0067 # molecular mass (kg/kmol)
        
        if self.tank_ox_material == "Aluminum 2024-T6":
            a = 23.4 * pow(10, -6)  # (1/K)
            self.rho_ox_tank = 2700  # kg/m^3
        elif self.tank_ox_material == "Aluminum 2014-T6":
            a = 20.8 * pow(10,-6) # (1/K)
            self.rho_ox_tank = 2800 # kg/m^3

        if self.tank_f_material == "Aluminum 2024-T6":
            a = 23.4 * pow(10, -6)  # (1/K)
            self.rho_f_tank = 2700  # kg/m^3
        elif self.tank_f_material == "Aluminum 7075-T6":
            a = 23.6 * pow(10, -6)  # (1/K)
            self.rho_f_tank = 2810  # kg/m^3
        elif self.tank_f_material == "Aluminum 2014-T6":
            a = 20.8 * pow(10,-6) # (1/K)
            self.rho_f_tank = 2800 # kg/m^3

        if self.fairing_material == "Aluminum 6061-T6":
            self.rho_fairing = 2700 # kg/m^3  aluminum 6061 density
            #self.fairing_Ftu = 310 # MPa
        elif self.fairing_material == "Aluminum 7075-T6":
            self.rho_fairing = 2810 # kg/m^3 al 7075-T6 density
        elif self.fairing_material == "Aluminum 2219-T852":
            self.rho_fairing = 2810 # kg/m^3 al 7075-T6 density

        # Parameters for actual tank volumes
        step_dia = 2*self.r  # diameter (m)
        self.circumference = 2*pi*self.r  # m Circumference

        # Propellant Mass Calculations
    
        self.residual_prop_perc = 0.02  # residual propellant (1)

        if (self.LV.Mission.input[0] == 'One') | (self.step_num != self.LV.num_steps):
            self.T_SL = self.LV.m_0[self.step_num-1]*g*self.TW
        elif self.LV.name == 'Latona-2':
            self.T_SL = (self.LV.m_0[self.step_num-1]-65)*g*self.TW
        elif (self.LV.name == 'Minerva-2') | (self.LV.name == 'Zephyr-2'):
            self.T_SL = (self.LV.m_0[self.step_num-1] + 30)*g*self.TW

        # if self.LV.name == 'Latona-1':
        #     self.T_SL = self.LV.m_0[self.step_num-1]*g*self.TW
        # elif self.LV.name == 'Latona-2':
        #     if self.step_num == 3:
        #         self.T_SL = (self.LV.m_0[self.step_num-1]-65)*g*self.TW
        #     else:
        #         self.T_SL = self.LV.m_0[self.step_num-1]*g*self.TW
        # elif self.LV.name == 'Minerva-1':
        #     self.T_SL = self.LV.m_0[self.step_num-1]*g*self.TW
        # elif self.LV.name == 'Minerva-2':
        #     if self.step_num == 2:
        #         self.T_SL = (self.LV.m_0[self.step_num-1] + 30)*g*self.TW
        #     else:
        #         self.T_SL = self.LV.m_0[self.step_num-1]*g*self.TW
        #print(self.LV.name + " has stage " + str(self.step_num) + " thrust of " + str(self.T_SL) )
        if self.parallel:
            self.multiplier = self.num_boosters
        else: 
            self.multiplier =1
        # T_SL using T/W and backward iteration -> VERIFIED
        # if self.parallel: # Iterate backward through all steps to get T_SL from the T/W of GROSS lift-off mass
        #     self.T_SL = 0
        #     for i in range(len(self.LV.m_0))[::-1]:
        #         self.T_SL += self.LV.m_0[i]*g*(self.TW) # sum the total weights of each step and multiply by the T/W to get the Thrust supplied by this step
        #         if self.step_num -1  == i: # break after summing the current step's weight * T/W to T/SL
        #             break
        #     self.multiplier = self.num_boosters # set multiplier to the number of boosters
        # else:
        #     self.T_SL = 0
        #     for i in range(len(self.LV.m_0))[::-1]: # iterate backward through all steps until the current step to get T_SL from the T/W of the STAGE lift-off mass
        #         self.T_SL += self.LV.m_0[i]*g*self.TW # sum the total weights of each step and multiply by the T/W to get the Thrust supplied by this step
        #         if (self.LV.name == 'Minerva-2') & (self.step_num == 1):
        #             self.T_SL -= (self.LV.m_0[2]-65)*g*self.TW
        #             pass
        #         if self.step_num-1 == i: # break after summing the current step's weight * T/W to T/SL
        #             break
        #     self.multiplier = 1 # set multiplier to 1 since this step has no boosters
        # self.T_SL = self.T_SL / self.multiplier # divide by multiplier to size each step/booster accordingly
        self.T_SL_engine = self.T_SL/self.num_engines
        #disp(T_SL)
        self.residual_prop = self.residual_prop_perc * self.LV.m_p[self.step_num-1]/self.multiplier  # residual propellant per step/booster(kg)
        m_dot = self.T_SL / (Isp *g)  # mass flow rate per step/booster (kg/s)
        self.startup_prop = m_dot*self.t_start  # propellant needed for start-up per step/booster (kg)
        self.m_prop_tot = (self.residual_prop + self.startup_prop + self.LV.m_p[self.step_num-1])/self.multiplier # # total propellant needed per step/booster (kg)
        ## Needed Tank Volumes
        if self.propulsion == 'Liquid':
            self.m_f_ideal = self.m_prop_tot* self.fuel_frac
            self.vol_f_ideal = self.m_f_ideal/self.rho_f  # ideal fuel volume (kg)
            self.ox_frac = 1/(1+FO) 
            self.m_ox_ideal = self.m_prop_tot * self.ox_frac
            self.vol_ox_ideal = self.m_ox_ideal/self.rho_ox  # ideal ox volume (kg)
            ## Shrinkage
            # Shrinkage Parameters
            shrinkage_fuel = 3*a*(rm_temp - fuel_temp)  # shrinkage factor for fuel tank (1)
            shrinkage_ox = 3*a*(rm_temp - ox_temp)  # shrinkage factor for ox tank (1)
            #disp(shrinkage_ox)
            ## Total Tank Volumes
            # Tank Sizing Parameters
            self.ullage_frac = 0.025  # ullage factor (how much gas is in the propellant tank) (1)
            vol_tot_needed_f = self.vol_f_ideal * (1 + self.ullage_frac + shrinkage_fuel)
            vol_tot_needed_ox = self.vol_ox_ideal * (1 + self.ullage_frac + shrinkage_ox)
            # Pressurant tank calculation
            
            print("p_press is " + str(self.p_press))
            print("p_tank is " + str(p_tank))
            print("vol_tot_needed_ox is " + str(vol_tot_needed_ox))
            print("vol_tot_needed_f is " + str(vol_tot_needed_f))
            vol_press = p_tank * (vol_tot_needed_f + vol_tot_needed_ox)/(self.p_press - p_tank) # Calculate required pressurant tank volume
            print("vol_press is " + str(vol_press))
            r_press = pow(3 * vol_press/(4*pi), 1/3) # m calculate radius of pressure tank (spherical)
            # if self.pressure_fed:
            #     if r_press > 2*self.r: # size pressure tank to 2 * radius of the liquid step
            #         self.p_press += 0.101325 # add 1 atmosphere in MPa
            #         self.sizeStep()
            #         return
            #else:
            if r_press > self.r/2: # size pressure tank to 1/2 radius of the liquid step
                self.p_press += 0.101325 # add 1 atmosphere in MPa
                self.sizeStep()
                return

            self.press_tank.append(r_press) # m append radius of pressure tank dome
            SA_press_tank = 4*pi*pow(r_press,2) # m^2 calculate SA of pressure tank
            self.press_tank.append(SA_press_tank) # append surface area of pressure tank to list
            self.press_tank.append(vol_press) # append volume of pressure tank to list
            self.m_press = vol_press * self.p_press * pow(10,6) / (R_uni * ox_temp) * Mol_mass

            print("Radius of pressure tank dome for step number " + str(self.step_num) + " is " + str(r_press))
            print("pressurant tank pressure is " + str(self.p_press))
            print("Volume of the pressure tank is " + str (vol_press))
            # calculate thickness of pressure tank
            ## Actual Tank Volumes Calculations
            ## Fuel
            self.dome_f.append(self.r/AR)  # m Lengths of fuel dome (m)
            #disp(self.dome_f[0])
            self.dome_f.append(pi*pow(self.r,2)*( 1 + 1/(2*e*pow(AR,2))*log((1+e)/(1-e))))  # m^2 SA of the fuel dome (m^2)
            dome_vol_f = 2/3*pi*pow(self.r,2)*self.dome_f[0]  # m^3 dome fuel tank volume (m^3)
            self.dome_f.append(dome_vol_f) # append dome vol
            cyl_vol_f = vol_tot_needed_f - 2*dome_vol_f  # m^3 cylinder fuel tank volume (m^3)
            self.cyl_f.append(cyl_vol_f/(pi*pow(self.r,2)))  # m Lengths of the fuel cylinder (m)
            #print("This iteration has Fuel Cylinder height " +str(self.cyl_f[0]))
            #print("This iteration has step radius " + str(self.r))
            # Iterate sizeStep, changing step diameter to avoid a negative length of tank cylinder
            if self.cyl_f[0] < 0:
                self.r=self.r-0.01
                self.sizeStep() 
                return
     
            self.cyl_f.append(self.cyl_f[0]*self.circumference)  # m^2 SA of fuel cylinder
            self.cyl_f.append(cyl_vol_f) # append vol m^3 of fuel cylinder
            #disp(self.self.cyl_f[0])
            ## Oxidizer
            vol_tot_needed_ox += vol_press # add pressurant tank inside oxidizer propellant tank -> changes dimensions required to house both the 
            # oxidizer and the pressurant tank <-- DONE
            h_dome_ox = self.r/AR # Lengths of the oxidizer dome (m)
            self.dome_ox.append(h_dome_ox)  # append height of ox dome
            #disp(self.dome_ox[0])
            SA_dome_ox = pi*pow(self.r,2)*( 1 + 1/(2*e*pow(AR,2))*log((1+e)/(1-e))) # SA of the oxidizer dome (m^2)
            self.dome_ox.append(SA_dome_ox) # append SA of ox dome
            dome_vol_ox = 2/3*pi*pow(self.r,2)*self.dome_ox[0]  # dome oxidizer tank volume (m^3)
            self.dome_ox.append(dome_vol_ox) # append vol of ox dome
            cyl_vol_ox = vol_tot_needed_ox - 2*dome_vol_ox  # cylinder fuel tank volume (m^3)
            h_cyl_ox = cyl_vol_ox/(pi*pow(self.r,2)) # length of oxidizer cylinder
            self.cyl_ox.append(h_cyl_ox)  # append to list the length of the ox cylinder
            #disp(self.cyl_ox[0])

            # Iterate sizeStep, changing step diameter to avoid a negative length of tank cylinder
            if self.cyl_ox[0] < 0:
                self.r=self.r-0.01
                self.sizeStep()
                return
            SA_cyl_ox = h_cyl_ox*self.circumference
            self.cyl_ox.append(SA_cyl_ox)  # m^2 SA of oxidizer cylinder
            self.cyl_ox.append(cyl_vol_ox) # m^3 of oxidizer cylinder

            self.intertank.append(self.r/2+self.dome_f[0]+self.dome_ox[0])  # m length of intertank
            self.intertank.append(self.intertank[0]*self.circumference)   # m^2 surface area of intertank

            ## LIQUID BI-PROPELLANT ENGINE SIZING
            Lth_ratio = 0.9  # chamber length / chamber diameter ratio
            alpha_conv = radians(30)  # convergent section angle (radians)
            L_star = 0.75  # characteristic length (m)
            
            if self.step_num == self.LV.num_steps:
                T_vac = self.T_SL 
            else:
                T_vac = self.T_SL/(1-(p_inf*Isp*g*self.epsilon)/(self.p_c*g*Isp_vac))  # total vacuum thrust ACTUAL
                #T_vac = self.T_SL/(1-(p_inf*Isp*g*self.epsilon)/(self.p_c*g*311))  # total vacuum thrust ON EXCEL


            T_vac_engine = T_vac/self.num_engines  # vacuum thrust of each engine
            A_t = Isp*g*T_vac_engine/(self.p_c*g*Isp_vac)  # throat area ACTUAL
            #A_t = Isp*g*T_vac_engine/(self.p_c*g*Isp_vac*2)  # throat area ON EXCEL

            d_t = sqrt(4*A_t/pi)  # throat diameter
            self.d_e = self.epsilon*d_t  # exit diameter
            self.L_n = (self.d_e-d_t)*0.8/((2*tan(radians(30))))  # length of nozzle for 80% length & 15 deg cone
            #print(self.L_n)
            L_c = pow((4*L_star * A_t * Lth_ratio)/pi, 1/3)  # chamber length
            d_c=  L_c/Lth_ratio 
            L_conv = (d_c - d_t)/(2*tan(alpha_conv)) 
            #self.L_eng_tot = L_n + L_c + L_conv 
            #disp(L_conv)

        elif self.propulsion == 'Solid':
            self.k_SRM = 0.135*0.855 # k_SRM-GrE = 0.855 * k_SRM-steel = 0.1154
            #print(self.k_SRM)
            # Ideal volumes
            self.dome_f.append(self.r/AR)  # m Lengths of fuel dome (m)
            #print(self.dome_f[0])
            self.dome_f.append(pi*pow(self.r,2)*( 1 + 1/(2*e*pow(AR,2))*log((1+e)/(1-e))))  # m^2 SA of the fuel dome (m^2)
            dome_vol_f = 2/3*pi*pow(self.r,2)*self.dome_f[0]  # m^3 dome fuel tank volume (m^3)
            self.dome_f.append(dome_vol_f) # append dome vol

            self.m_prop_ideal = self.m_prop_tot
            self.vol_prop_ideal = self.m_prop_ideal / self.rho_prop

            # NO Shrinkage
            # NO ULLAGE
            burn_pattern = 0.12 # assume 12% of casing volume is empty space allocated for initial burn pattern
            vol_tot_needed_prop_casing = (self.vol_prop_ideal - 2*self.dome_f[2])* (1 + burn_pattern) # m^3 total volume needed for casing

            vol_tot_casing = vol_tot_needed_prop_casing # m^3 volume per SRM casing
            
            h_srm_casing = vol_tot_casing / (pi * pow(self.r,2)) # m, height of srm casing derived from V = pi * r^2 * h
            self.srm_casing.append(h_srm_casing) # append height to srm_casing list
            SA_srm_casing = h_srm_casing * pi * step_dia # m^2 SA of srm_casing
            self.srm_casing.append(SA_srm_casing) # append SA of srm_casing to list
            self.srm_casing.append(vol_tot_casing) # append vol of srm_sasing to list
            # ideal gas lP_1 * V_1 = P_2 * (V_1 + V+2) -> solve for V_1 = vol_press
            vol_press = p_tank * (vol_tot_casing)/(self.p_press - p_tank) # Calculate required pressurant tank volume
            
            # to fill the srm casing at the tank pressure, p_tank, with a pressurant pressure of self.p_press (20 MPa currently
            r_press = pow(3 * vol_press/(4*pi), 1/3) # m calculate radius of pressure tank (spherical)
            if self.pressure_fed:
                if r_press > 2*self.r: # size pressure tank to 2 * radius of the liquid step
                    self.p_press += 0.101325 # add 1 atmosphere in MPa
                    self.sizeStep()
                    return
            else:
                if r_press > self.r: # size pressure tank to radius of the liquid step
                    self.p_press += 0.101325 # add 1 atmosphere in MPa
                    self.sizeStep()
                    return
            #print("Mass of SRM Casing is " + str(self.m_prop_tot/self.multiplier))
            #print("Volume of SRM Casing is" + str(vol_tot_casing))
            self.press_tank.append(r_press) # m append length of pressure tank dome
            SA_press_tank = 4 * pi * pow(r_press,2) # m^2 calculate SA of pressure tank
            self.press_tank.append(SA_press_tank) # append surface area of pressure tank to list
            self.press_tank.append(vol_press) # append volume of pressure tank to list
            self.m_press = vol_press * self.p_press * pow(10,6) / (R_uni * rm_temp) * Mol_mass
            print("Radius of pressure tank dome is " + str(r_press))
            print("Pressurant tank pressure is " + str(self.p_press))
            print("Volume of the pressure tank is " + str(vol_press))
            ## Solid NOZZLE SIZING
            Lth_ratio = 0.9  # chamber length / chamber diameter ratio
            alpha_conv = radians(30)  # convergent section angle (radians)
            L_star = 0.75  # characteristic length (m)
            
            if self.step_num == self.LV.num_steps:
                T_vac = self.T_SL 
            else:
                T_vac = self.T_SL/(1-(p_inf*Isp*g*self.epsilon)/(self.p_c*g*Isp_vac))  # total vacuum thrust ACTUAL
                #T_vac = self.T_SL/(1-(p_inf*Isp*g*self.epsilon)/(self.p_c*g*311))  # total vacuum thrust ON EXCEL


            T_vac_engine = T_vac/self.num_engines  # vacuum thrust of each engine
            self.L_n = self.r
            self.d_e = self.r
        ## Intertank and Thrust Structure Sizing
        # The following variables denote length
        
        self.getSkirtsnPL()
        self.getTotLength()
      
    def getTotLength(self):
    # get total length of step
        if self.step_num == self.LV.num_steps:
            if self.propulsion == 'Liquid':
                self.total_length = self.aft_skirt[0] + self.fwd_skirt[0] + self.intertank[0] + self.cyl_f[0] + self. cyl_ox[0] 
            elif self.propulsion == 'Solid':
                self.total_length = self.aft_skirt[0] + self.fwd_skirt[0] + self.srm_casing[0]
        else:
            if self.propulsion == 'Liquid':
                self.total_length = self.aft_skirt[0] + self.intertank[0] + self.cyl_f[0] + self. cyl_ox[0] 
            elif self.propulsion == 'Solid':
                self.total_length = self.aft_skirt[0] + self.srm_casing[0]

    def getSkirtsnPL(self):
        if self.step_num == self.LV.num_steps:
            self.fwd_skirt.append(2/3*self.r + self.dome_ox[0])  # m Length of fwd skirt
            self.fwd_skirt.append(self.fwd_skirt[0]*self.circumference)  # m^2 SA of fwd skirt
            self.aft_skirt.append(self.fwd_skirt[0])  # m Length of aft skirt
            self.aft_skirt.append(self.aft_skirt[0]*self.circumference)  # m^2 SA of aft skirt
            # NOTE: Change T_struct based on nozzle length
            self.T_struct = self.aft_skirt[0]/2  # length of thrust structure
            

            # Payload VOLUME --> VERIFY WITH SUNNY
            #cubeDis_length = 0.13  # m sidelength of cubeSat dispenser
            if self.LV.name == "Latona-1":
                h_fairing_cyl = 0.05 + 0.414 # from solidworks layout -> 4 3U's house 12 cubesats -> 30 cm height; add 4 cm for thickness and lip
            elif self.LV.name == "Latona-2":
                h_fairing_cyl = 0.01 + 0.45 + 0.35 # from solidworks layout -> 5 3U's house 15 cubesats, 4 1U's and 1 2U house the other 6 -> 50 cm height; add 4 cm for thickness and lip
            elif (self.LV.name == "Minerva-1") | (self.LV.name == "Zephyr-1"):
                h_fairing_cyl = 0.01 + 0.45 # from solidwroks layout single platform width -> 2 cm for thickness and lip
            elif (self.LV.name == "Minerva-2") | (self.LV.name == "Zephyr-2"):
                h_fairing_cyl = 0.01 + 0.45 # from solidworks layout, 2 platforms -> 20 cm height, add 4 cm for thickness and lip

            if self.fairing_shape == "Cone":
                h_fairing_cone = self.r*2
                self.fairing.append([h_fairing_cyl, h_fairing_cone]) # m  [fairing cylinder height, fairing cone height]
                l = pow(pow(h_fairing_cone, 2) + pow(self.r,2), 1/2) # l parameter for calculating cone surface area
                self.fairing.append([2*(self.r)*pi*h_fairing_cyl, pi*self.r*l])  # m^2 fairing surface area
                # larger payload radius (0.275 m) than step radius (0.23 m)
            #print("fairing cylinder height is " + str(self.fairing[0][0]))
            #print("fairing cone height is " + str(self.fairing[0][1]))
            print("fairing cylinder SA is " + str(self.fairing[1][0]))
            print("fairing cone SA is " + str(self.fairing[1][1]))
            #print("payload mass is " + str(self.m_PL))
        elif self.step_num != self.LV.num_steps:
            self.aft_skirt.append(2*self.r)  # m length of the aft skirt
            self.aft_skirt.append(self.aft_skirt[0]*self.circumference) 
            # NOTE: Change T_struct based on nozzle length
            self.T_struct = self.aft_skirt[0]/2  # m length of thrust structure

            if self.parallel:
                if self.fairing_shape == "Cone":
                    h_fairing_cone = self.r*2
                    self.fairing.append(h_fairing_cone)
                    l = pow(pow(h_fairing_cone, 2) + pow(self.r,2), 1/2)
                    self.fairing.append(pi*self.r*l)  # m^2 fairing surface area for cone
                #print("Booster nose cone height is " + str(self.fairing[0]))

    def setInsulation(self):
        self.SA_rho_insulation = {'Oxidizer': 0, 'Fuel': 0}
        if self.insulation_dict['Fuel'] == 'Rubber':
            self.SA_rho_insulation['Fuel'] = 1.123 # kg/m^2 LOX insulation 
        if self.insulation_dict['Oxidizer'] == 'Rubber':
            self.SA_rho_insulation['Oxidizer'] = 1.123 # kg/m^2 LOX insulation