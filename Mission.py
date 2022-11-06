import math
import numpy as np
from math import sqrt, pow, pi, cos, sin, tan, acos, asin, atan, radians, degrees, exp, log

class Mission:
    """Declares the Mission parameters inside a Mission object. Uses function getTrajReqs to find the delta-v trajectory requirements needed to size the LV"""
    
    # CONSTANTS OF EARTH
    g = 0.00980065  # gravity (km/s^2)
    mu_E = 398600  # gravitational parameter of earth (km^3/s)
    r_E = 6378  # radius of earth (km)
    v_equator = 0.4651 # equatorial velocity in (km/s)
    dV_reqs_names = ['dv needed', 'dv design', 'dv plane change', 'dv gravity loss', 'dv drag loss', 'dv apo kick', 'dv maneuvers']
    # Delta-V Trajectory Requirements Variable
    dV_reqs = '\'dv_reqs\' is not yet initialized, call the set method \'set_dV_reqs\''
    # CONSTRUCTOR
    def __init__(self, Mission_type, recovery, losses_gravity, drag_loss, launch_site):
        #Mission_type Construct an instance of this class
        #   Detailed explanation goes here
        self.input = [Mission_type, recovery, losses_gravity, drag_loss, launch_site]
        #self.dV_reqs_names = ['dv needed', 'dv design', 'dv plane change', 'dv gravity loss', 'dv drag loss', 'dv apo kick', 'dv maneuvers']

    # Delta-V Design Function
    def set_dV_reqs(self):
        #print(self.input)
        if self.input[0] == 'One':
            self.payload = 30
            delta_plane = 0
            inc = radians(60)
            h_a = 500 # apoapsis altitude (km)
            h_p = 200 # periapsis altitude (km)
        elif self.input[0] =='Two':
            self.payload = 95
            delta_plane = radians(10) # plane change of 10 degrees (rad)
            inc = radians(98) # inclination (radians)
            h_a = 550 # apoapsis altitude (km)
            h_p = 200 # periapsis altitude (km)
        elif self.input[0] =='Three':
            self.payload = 95
            delta_plane = radians(10) # plane change of 10 degrees (rad)
            inc = radians(98) # inclination (radians)
            h_a = 550 # apoapsis altitude (km)
            h_p = 200 # periapsis altitude (km)
        elif self.input[0] == 'LEAP':
            self.payload = 1
            delta_plane = 0
            inc = radians(0)
            h_a = 15.24
            h_p = 15.24 # periapsis altitude (km)

        # h_p = 200 # periapsis altitude (km)
        if self.input[4] == 'Kodiak':
            self.lat = radians(57.8324683) # latitude (radians)

        elif self.input[4] == 'Vandenberg':
            self.lat = radians(34.7331518097343)
   
        elif self.input[4] == 'KSC':
            self.lat = radians(28.5226326595524)
        
   
        ### Orbital Calculations
        #  Orbital velocities
        r_p = self.r_E + h_p # periapsis radius (km)
        r_a = self.r_E + h_a # apoapsis radius (km)
        a = (r_p + r_a)/2 # semi-major axis (km)
        #print(r_p, r_a, a)

        v_p = sqrt(2*self.mu_E*(1/r_p - 1/(2*a))) # periapsis velocity (km/s)
        #print(v_p)
        v_a = sqrt(2*self.mu_E*(1/r_a - 1/(2*a))) # apoapsis velocity (km/s)
        #print(v_a)
        v_c = sqrt(self.mu_E/r_a) # circular velocity at h_a (km/s)
        #print(v_c)
        v_LS = self.v_equator * cos(self.lat) # launch site velocity (km/s)
        #print(v_p, v_a, v_c, v_LS)

    

        # Angles
        if inc > pi/2:
            aux = pi - inc # launch window auxiliary angle (rad)
        elif inc < pi/2:
            aux = inc # launch window auxiliary angle (rad)
        flt_path = asin(cos(aux)/cos(self.lat)) # flight path angle (rad)
        #print(aux, flt_path)

        # Azimuth Angle and Burnout Velocities
        if inc > pi/2:
            azimuth = pi + flt_path # azimuth angle (rad)
            v_BO_S = -v_p*cos(flt_path)*cos(azimuth) # South burnout velocity (km/s)
            v_BO_E = -v_p*cos(flt_path)*sin(azimuth) # East burnout velocity (km/s)
            v_BO_Z = v_p*sin(flt_path) # Zenith burnouth velocity (km/s)

        elif inc < pi/2:
            azimuth = flt_path
            v_BO_S = -v_p*cos(flt_path)*cos(azimuth) # South burnout velocity (km/s)
            v_BO_E = v_p*cos(flt_path)*sin(azimuth) # East burnout velocity (km/s)
            v_BO_Z = v_p*sin(flt_path) # Zenith burnouth velocity (km/s)

        #print(v_BO_S, v_BO_E, v_BO_Z)
        v_N_S = v_BO_S # South needed delta-v (km/s)
        v_N_E = v_BO_E - v_LS # East needed delta-v (km/s)
        v_N_Z = v_BO_Z # Zenith needed delta-v (km/2)
        #print(v_N_S, v_N_E, v_N_Z)

        dv_maneuvers = 0
        if self.input[1]:
            dv_maneuvers += 0.343 # landing dv burn (km/s)
        else: dv_maneuvers +=0

        dv_N = sqrt(pow(v_N_S,2) + pow(v_N_E,2) + pow(v_N_Z,2)) # total delta-v needed (km/s)
        grav_loss = 0
        if (type(self.input[2]) == str) & (self.input[2] == '80% gravity loss'):
            grav_loss = 0.8*sqrt(2*self.mu_E*h_p/((h_p+self.r_E)*self.r_E)) # 80# gravity loss eqn (km/s)
        elif (type(self.input[2]) == int) | (type(self.input[2]) == float):
            grav_loss = self.input[1]
        apo_kick = v_c - v_a # apoapsis kick burn (km/s)

        if self.input[0] == 'One':
            dv_design = dv_N + grav_loss + self.input[3] + apo_kick + dv_maneuvers # delta-v design is the total delta v required (km/s)
            #print(dv_N, grav_loss, drag_loss, apo_kick, dv_maneuvers, dv_design)
            self.dV_reqs = [dv_N, dv_design, 0, grav_loss, self.input[3], apo_kick, dv_maneuvers]
        elif self.input[0] =='Two':
            dv_plane = 2*v_a*sin(delta_plane/2) # delta-v needed for plane-change (km/s)
            dv_design = dv_N + grav_loss + self.input[3] + apo_kick + dv_maneuvers + dv_plane # delta-v design is the total delta v required (km/s)
            #print(dv_N, grav_loss, drag_loss, apo_kick, dv_maneuvers, dv_plane, dv_design)
            self.dV_reqs = [dv_N, dv_design, dv_plane, grav_loss, self.input[3], apo_kick, dv_maneuvers]
        elif self.input[0] =='Three':
            dv_plane = 2*v_a*sin(delta_plane/2) # delta-v needed for plane-change (km/s)
            dv_design = dv_N + grav_loss + self.input[3] + apo_kick + dv_maneuvers + dv_plane # delta-v design is the total delta v required (km/s)
            #print(dv_N, grav_loss, drag_loss, apo_kick, dv_maneuvers, dv_plane, dv_design)
            self.dV_reqs = [dv_N, dv_design, dv_plane, grav_loss, self.input[3], apo_kick, dv_maneuvers]
    
    def print(self):
        print('Here are the mission requirements for Mission ' + self.input[0] + ', where the recovery is ' + str(self.input[1]) + ' and the launch site is ' + str(self.input[4]))
        print("The delta-v needed is " + str(self.dV_reqs[0]) + " km/s")
        print("The delta-v design is " + str(self.dV_reqs[1]) + " km/s")
        print("The delta-v plane change is " + str(self.dV_reqs[2]) + " km/s")
        print("The delta-v gravity loss is " + str(self.dV_reqs[3]) + " km/s")
        print("The delta-v drag loss is " + str(self.dV_reqs[4]) + " km/s")
        print("The delta-v apo-kick is " + str(self.dV_reqs[5]) + " km/s")
        print("The delta-v maneuvers (entry and landing burn) is " + str(self.dV_reqs[6]) + " km/s")
        print()

    # def GoalSeek(self, fun,goal,x0,fTol=0.0001,MaxIter=1000):
    # # Goal Seek function of Excel
    # #   via use of Line Search and Bisection Methods

    # # Inputs
    # #   fun     : Function to be evaluated
    # #   goal    : Expected result/output
    # #   x0      : Initial estimate/Starting point

    # # Initial check
    #     if fun(x0)==goal:
    #         print('Exact solution found')
    #         return x0

    #     # Line Search Method
    #     step_sizes=np.logspace(-1,4,6)
    #     scopes=np.logspace(1,5,5)

    #     vFun=np.vectorize(fun)

    #     for scope in scopes:
    #         break_nested=False
    #         for step_size in step_sizes:

    #             cApos=np.linspace(x0,x0+step_size*scope,int(scope))
    #             cAneg=np.linspace(x0,x0-step_size*scope,int(scope))

    #             cA=np.concatenate((cAneg[::-1],cApos[1:]),axis=0)

    #             fA=vFun(cA)-goal

    #             if np.any(np.diff(np.sign(fA))):

    #                 index_lb=np.nonzero(np.diff(np.sign(fA)))

    #                 if len(index_lb[0])==1:

    #                     index_ub=index_lb+np.array([1])

    #                     x_lb=np.asscalar(np.array(cA)[index_lb][0])
    #                     x_ub=np.asscalar(np.array(cA)[index_ub][0])
    #                     break_nested=True
    #                     break
    #                 else: # Two or more roots possible

    #                     index_ub=index_lb+np.array([1])

    #                     print('Other solution possible at around, x0 = ', np.array(cA)[index_lb[0][1]])

    #                     x_lb=np.asscalar(np.array(cA)[index_lb[0][0]])
    #                     x_ub=np.asscalar(np.array(cA)[index_ub[0][0]])
    #                     break_nested=True
    #                     break

    #         if break_nested:
    #             break
    #     if not x_lb or not x_ub:
    #         print('No Solution Found')
    #         return

    #     # Bisection Method
    #     iter_num=0
    #     error=10

    #     while iter_num<MaxIter and fTol<error:
        
    #         x_m=(x_lb+x_ub)/2
    #         f_m=fun(x_m)-goal

    #         error=abs(f_m)

    #         if (fun(x_lb)-goal)*(f_m)<0:
    #             x_ub=x_m
    #         elif (fun(x_ub)-goal)*(f_m)<0:
    #             x_lb=x_m
    #         elif f_m==0:
    #             print('Exact spolution found')
    #             return x_m
    #         else:
    #             print('Failure in Bisection Method')
        
    #         iter_num+=1

    #     return x_m



