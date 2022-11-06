# SCRIPT SUMMARY:
# This is the main script used to generate mass budgets and size launch vehicles to complete a mission for a given altitude, inclination, and payload.
# First the launch vehicle's propellant mass budgets are estimated by calling a MATLAB script in the function "generateLVMassEstimates(missions)", which takes a list of mission objects as the parameter.
# Mission objects detail the orbital requirements and the required delta-v to complete the mission.
# After the propellant mass budgets are calculated, the "LaunchVehicle" objects are initialized, which are the parents of the "Step" objects.
# The "Step" objects initialize parameters such as propellant choice, gas used to pressurize tanks, materials, dome shape, and number of engines. The geometry of the step's body and tank components are then sized to accomodate the required volume of propellant.
# After each "Step" is sized, the interstages of the "LaunchVehicle" are sized and the moments of inertia of each launch vehicle component are tabulated into a table.
# A MATLAB script is called in "runTrajectory()" to run trajectory simulations of the "LaunchVehicle" objects using a gravity-turn and the optimal trajectory is selected based on quickest time to orbit and largest amount of propellant left for reserves.
# The conditions at the point of maximum dynamic pressure (Max-Q) for the optimal trajectory is tabulated, such as: angle of attack, vehicle velocity, wind speed, left-over propellant mass.
# SCRIPT DOES NOT YET AUTOMATE:::
# Using the Max-Q data, Axial Load, Shear Load, and Bending Moment diagrams are made for the "LaunchVehicle" at the Max-Q condition, and Ground Wind Load Condition in Excel.
# From there, the stresses through the "LaunchVehicle" are tabulated in Excel, and the thicknesses of body components are selected to prevent stresses at the 
# Ultimate Compressive Strength of the material by at least a margin of 50%. This depends on which components of the "LaunchVehicle" are pressurized and which ones are not.
# The script also exports a PPT presentation slide of a to-scale, side-view profile of the "LaunchVehicle". <-- Work still must be done to position PPT shapes correctly

import math
from math import sqrt, pow, pi, cos, sin, acos, asin, atan, radians, degrees
from Mission import Mission
from LaunchVehicle import LaunchVehicle
from Step import Step
from pptx import Presentation
import matlab.engine
import pandas as pd

def generateLVMassEstimates(missions):
    # This function loops through the passed in argument 'missions', which is a
    #  list of Mission objects, and for each Mission object i, a pandas data-
    # frame is contructed  
    # missions variable is mission.dV_reqs
    print("Generating LV mass estimates through MATLAB...")
    for i in missions:
        df = pd.DataFrame(columns=[i.dV_reqs_names[0], i.dV_reqs_names[1], i.dV_reqs_names[2], i.dV_reqs_names[3],
                            i.dV_reqs_names[4], i.dV_reqs_names[5], i.dV_reqs_names[6], 'Latitude', 'Payload'], index=range(1))
        for j in range(len(i.dV_reqs_names)):
            df[i.dV_reqs_names[j]][0] = i.dV_reqs[j]*1000
        df.iloc[0, len(i.dV_reqs_names) + 1] = i.payload
        df.iloc[0, len(i.dV_reqs_names)] = i.lat
        df.to_csv('LVMasses\\' + i.input[0] + 'dVParameters.csv', index=False)
        pass
        # Generate csv file for each mission in the list 'missions' that contains all pertinent parameters required by MATLAB script
    eng = matlab.engine.start_matlab()
    #eng.Mass_Estimates_PythonLinked(nargout=0)
    eng.Mass_Estimates_Zephyr_PythonLinked(nargout=0)
    print("Mass estimates generated.")

def runTrajectory():
    print("Running MATLAB Trajectory...")
    eng = matlab.engine.start_matlab()
    #eng.Trajectory_Run_TL_0324(nargout=0) # For Latona and Minerva
    eng.Trajectory_Run_TL_0407(nargout=0) # For Zephyr
    print("MATLAB Trajectory Complete.")

material = ('Aluminum 6061-T6', 'Rubber', 'Aluminum 2024-T6', 'Aluminum 2014-T6', 'Aluminum 7075-T6', 'Aluminum 2219-T87', 'Aluminum 2219-T852') # materials
grav_est = ('80% gravity loss', 0.5, 1, 1.5, 2) # gravity estimates (km/s)
drag_est = (0.2, 0.4, 0.6, 0.8, 1.0) # drag estimates (km/s) 

# insulation dictionaries for each step in Latona 1
no_insul = {'Fuel': None, 'Oxidizer': None} # no insulation for oxidizer or fuel tank
ox_insul = {'Fuel': None, 'Oxidizer': material[1]} # insulation for only oxidizer
all_insul = {'Fuel': material[1], 'Oxidizer': material[1]} # insulation for oxidizer and fuel tanks

launch_site = ('Kodiak', 'KSC', 'Vandenberg') # launch site locations

missions = []
One = Mission('One', False, grav_est[0], drag_est[0], launch_site[1])
missions.append(One)
Two = Mission('Two', False, grav_est[0], drag_est[0], launch_site[2])
missions.append(Two)
One.set_dV_reqs()
Two.set_dV_reqs()
generateLVMassEstimates(missions)

engine = ('Raptor', 'Merlin', 'SRBNozzle') # engine type to be used for scaling
propellant = ('Kerolox', 'Methalox', 'AP-Al-HTPB') # propellant type
pressurant = ('Helium', 'Nitrogen')
dome_shape = ('Circular', 'Elliptical-2', 'Elliptical-sqr2') # propellant tank dome shape
fairing_shape = ('Cone', 'Ogive') # payload fairing shape
loads_conditions = ('Ground Wind-Loads Condition', 'Max-Q Condition')
#print(engine.index('Merlin'))
LVSteps = []
LaunchVehicles = []

Latona1 = LaunchVehicle('Latona-1', 1.4, material[0], 2, [265, 380], [0.11, 0.11], 30, One)
LaunchVehicles.append(Latona1)
                  #(LV, r, step_num, TW, engine, dome shape, propellant, pressurant, insulation_dict, tank material[ox,f], number of engines, num gimballed engines, t_start, fairing_material, fairing_shape, num_boosters)
L1_step1 = Step(Latona1, 0.275, 1, 1.4, engine[1], dome_shape[1], propellant[2], pressurant[0], no_insul, [material[3], material[0]], 1, 1, 1, None, None, 0)
L1_step2 = Step(Latona1, 0.275, 2, 0.9, engine[0], dome_shape[1], propellant[1], pressurant[0], all_insul, [material[3],material[3]], 1, 1, 1, material[6], fairing_shape[0], 0)

Latona_steps = [L1_step1, L1_step2]
LVSteps.append(Latona_steps)

Latona2 = LaunchVehicle('Latona-2', 2.0, material[0], 3, [265, 265, 380], [0.11, 0.11, 0.11], 95, Two)
LaunchVehicles.append(Latona2)
L2_step1 = Step(Latona2, 0.275, 1, 1, engine[1], dome_shape[1], propellant[2], pressurant[0], no_insul, [material[3], material[0]], 1, 0, 1, material[6], fairing_shape[0], 4)     # boosters
L2_step2 = Step(Latona2, 0.275, 2, 1.4, engine[1], dome_shape[1], propellant[2], pressurant[0], no_insul, [material[3], material[0]], 1, 1, 1, None, None, 0)                       # main stage
L2_step3 = Step(Latona2, 0.275, 3, 0.9, engine[0], dome_shape[1], propellant[1], pressurant[0], all_insul, [material[3],material[3]], 1, 1, 1, material[6], fairing_shape[0], 0)   # upper stage

Latona_steps = [L2_step1, L2_step2, L2_step3]
LVSteps.append(Latona_steps)

Minerva1 = LaunchVehicle('Minerva-1', 1.4, material[4], 2, [296.1, 359.1], [0.11, 0.11], 30, One)
LaunchVehicles.append(Minerva1)
#                    # (r, laststep, TW, engine,  dome shape,propellant, body material, tank material, number of engines, t_start, fairing_material, fairing_shape) 
M1_step1 = Step(Minerva1, 0.45, 1, 1.4, engine[1], dome_shape[1], propellant[0], pressurant[0], ox_insul, [material[3], material[4]], 1, 1, 1, None, None, 0)
M1_step2 = Step(Minerva1, 0.45, 2, 1.05, engine[1], dome_shape[1], propellant[0], pressurant[0], ox_insul, [material[3], material[4]], 1, 1, 1, material[6], fairing_shape[0], 0)
Minerva_steps = [M1_step1, M1_step2]
LVSteps.append(Minerva_steps)

Minerva2 = LaunchVehicle('Minerva-2', 1.2, material[4], 3, [296.1, 359.1, 359.1], [0.11, 0.11, 0.11], 95, Two)
LaunchVehicles.append(Minerva2)
#                    # (r, laststep, TW, engine,  dome shape,propellant, body material, tank material, number of engines, t_start, fairing_material, fairing_shape) 
M2_step1 = Step(Minerva2, 0.45, 1, 1.4, engine[1], dome_shape[1], propellant[0], pressurant[0], ox_insul, [material[3], material[4]], 1, 1, 1, None, None, 0)
M2_step2 = Step(Minerva2, 0.45, 2, 1.05, engine[1], dome_shape[1], propellant[0], pressurant[0], ox_insul, [material[3], material[4]], 1, 1, 1, None, None, 0)
M2_step3 = Step(Minerva2, 0.39, 3, 0.9, engine[1], dome_shape[1], propellant[0], pressurant[0], ox_insul, [material[3], material[4]], 1, 1, 1, material[6], fairing_shape[0], 0)
Minerva_steps = [M2_step1, M2_step2, M2_step3]
LVSteps.append(Minerva_steps)

Zephyr1 = LaunchVehicle('Zephyr-1', 1.4, material[0], 2, [330, 380], [0.11, 0.11], 30, One)
LaunchVehicles.append(Zephyr1)
#                    # (r, laststep, TW, engine,  dome shape,propellant, body material, tank material, number of engines, t_start, fairing_material, fairing_shape, pressure fed) 
Z1_step1 = Step(Zephyr1, 0.45, 1, 1.4, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, None, None, 0, True)
Z1_step2 = Step(Zephyr1, 0.45, 2, 1.05, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, material[6], fairing_shape[0], 0, True)
Zephyr1_steps = [Z1_step1, Z1_step2]
LVSteps.append(Zephyr1_steps)

Zephyr2 = LaunchVehicle('Zephyr-2', 1.3, material[0], 3, [330, 380, 380], [0.11, 0.11, 0.11], 95, Two)
LaunchVehicles.append(Zephyr2)
#                    # (r, laststep, TW, engine,  dome shape,propellant, body material, tank material, number of engines, t_start, fairing_material, fairing_shape) 
Z2_step1 = Step(Zephyr2, 0.45, 1, 1.4, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, None, None, 0, True)
Z2_step2 = Step(Zephyr2, 0.45, 2, 1.05, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, None, None, 0, True)
Z2_step3 = Step(Zephyr2, 0.39, 3, 0.9, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, material[6], fairing_shape[0], 0, True)
Zephyr2_steps = [Z2_step1, Z2_step2, Z2_step3]
LVSteps.append(Zephyr2_steps)

#print(LaunchVehicles)
for i in range(len(LaunchVehicles)):
    LV = LaunchVehicles[i]
    LV.initMassEstimates()
    LV.initSteps(LVSteps[i])
    for j in LV.listOfSteps:
        j.sizeStep()
    LV.initInterstages()
    LV.massMoments(loads_conditions[0])
    LV.addSlide()
    LV.generateTrajReqs()
runTrajectory() # function runs trajectory for launch vehicles specified inside function
for i in range(len(LaunchVehicles)):
    LV = LaunchVehicles[i]
    LV.massMoments(loads_conditions[1]) # wind loads
print('The program is finished!')