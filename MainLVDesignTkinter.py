# NOTES (Indented are completed):
    # Must still add gap between the top of the solid propellant casing and the next step for the pressure tank <-- DONE
    # Must still re-adjust total volume needed in the cryogenic tank because the pressure tank must be added and suspended inside of it <-- DONE
    # Add Nose Cone to boosters < -- DONE
    # ^ Calculate needed amount of pressurant to fill itself and the propellant tank at 30-60 psi <-- DONE
    # Re-calculate the J0's for the propellant tanks and insulation to account for dome <-- DONE
    # Recalculate mass estimating relationships for tanks, srm_casing <-- DONE
    # Distance calcs: place avionics in the forward skirt <-- DONE
    # Calculate payload fairing sizing by multiplying the # of CubeSats by its volume (0.001 m^3), multiply by 1.5 for effective volume to include <-- DONE
    # the attach fitting in the volume and other inefficiencies in volume allocation
    # Fix Payload fairing distances CG calculation to account for cylinder AND cone portions of the payload fairing
# Try to fix Scale Profile Drawing Views and upload them to the OneDrive to help Sunny and Dario dimension the launch vehicle properly
    # Change tank material
    # Change payload fairing from cone to spherically blended tangent ogive
    # Pick and justify a material - Al 2024, 5083, or 7039 for cryogenic tanks, 6061 for body
    # Aluminum 2024 - most popular in aerospace engineering -> allows small amounts of cold
    # deformation and increases yield strength. Used primarily in sheet forms for fuselage
    # and wings and has an Ftu of 470 MPa.
    # Stage 1 Minerva - Aluminum 7075 T6 - High strength and fatigue strength, easy machinability, resistant to corrosion https://www.aerospacemanufacturinganddesign.com/article/aluminum-alloys-for-aerospace/  https://www.engineeringclicks.com/7075-t6-aluminium/
    # Cryogenics - Aluminum 2014 - good for cryogenics, good machining, welding, low cost, APPLICATION: Cryo Tank Walls https://www.clintonaluminum.com/aluminum-in-cryogenic-applications/
    # All stages? - Aluminum 6061 T6 for body - higher (than 2014) strength, better welding and machining than 7075, low cost, https://parts-badger.com/6061-vs-7075-aluminum/
    # Fairing - Aluminum 2219-T852 - Has the highest operating temperature among aluminum alloys, lightweight, good weldability, strong, required heat treat for welds to resist corrosion https://www.aerospacemanufacturinganddesign.com/article/aluminum-alloys-for-aerospace/
    # Check Latona 1's propellant mass on LV Mass Moments vs what is given by Matlab-generated .csv, LVMasses. Check this for all LVs
    # Decrease Wiring Mass by a factor of 10?
    # Run Chase's LV Families MATLAB script using python; rewrite script to export csv's to 'C:\Dev\Visual Studios Solutions\Python\LaunchVehicleDesign\LVMasses'
    # Download Chase's Ground and Wind Loads Excel files and edit the sheets by rewriting the LVMassMoments table with the python generated LVMassMoments(dataframe)
    # ^ IMPORTANT: Must first reorder the stages to go Payload > Stage 3 > Stage 2 < Stage 1, instead of Payload > Stage 1, Stage 2, Stage 3
    # For automating trajectory run: Tyler must manually change the values of the variables for stage masses -> export those values in csv format and when refer
    # to that csv when initializing those variables in the Matlab Script
    # Get rid of thrust structure for SRMs, add tank domes on top and bottom of the srm casing

import math
from math import sqrt, pow, pi, cos, sin, acos, asin, atan, radians, degrees
from Mission import Mission
from LaunchVehicle import LaunchVehicle
from Step import Step
from pptx import Presentation
from tkinter import *
import matlab.engine
import pandas as pd

def generateLVMassEstimates(missions):
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
    # Edit Matlab file 'Mass_Estimates' to pull parameters from the csv files where it is currently being manually written.
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
grav_est = ('80% gravity loss', 0.5, 1, 1.5, 2) # gravity estimates (km)
drag_est = (0.2, 0.4, 0.6, 0.8, 1.0) # drag estimates (km) 

# insulation dictionaries for each step in Latona 1
no_insul = {'Fuel': None, 'Oxidizer': None} # no insulation for oxidizer or fuel tank
ox_insul = {'Fuel': None, 'Oxidizer': material[1]} # insulation for only oxidizer
all_insul = {'Fuel': material[1], 'Oxidizer': material[1]} # insulation for oxidizer and fuel tanks

launch_site = ('Kodiak', 'KSC', 'Vandenberg') # launch site locations

# Initialize Missions
missions = []
# One = Mission('One', False, grav_est[0], drag_est[0], launch_site[1])
# missions.append(One)
# Two = Mission('Two', False, grav_est[0], drag_est[0], launch_site[2])
# missions.append(Two)
#Three = Mission ('Three', False, grav_est[0], drag_est[0], launch_site[0])
# One.set_dV_reqs()
# Two.set_dV_reqs()
# Three.set_dV_reqs()
# Three.print()
#One.print()
#Two.print()

# Generate Mass Estimate .csv's for LVs in LVMasses folder
# generateLVMassEstimates(missions)

engine = ('Raptor', 'Merlin', 'SRBNozzle') # engine type to be used for scaling
propellant = ('Kerolox', 'Methalox', 'AP-Al-HTPB') # propellant type
pressurant = ('Helium', 'Nitrogen')
dome_shape = ('Circular', 'Elliptical-2', 'Elliptical-sqr2') # propellant tank dome shape
fairing_shape = ('Cone', 'Ogive') # payload fairing shape
loads_conditions = ('Ground Wind-Loads Condition', 'Max-Q Condition')
#print(engine.index('Merlin'))
LVSteps = []
LaunchVehicles = []

# Zephyr1 = LaunchVehicle('Zephyr-1', 1.4, material[0], 2, [330, 380], [0.11, 0.11], 30, One)
# LaunchVehicles.append(Zephyr1)
# #                    # (r, laststep, TW, engine,  dome shape,propellant, body material, tank material, number of engines, t_start, fairing_material, fairing_shape) 
# Z1_step1 = Step(Zephyr1, 0.45, 1, 1.4, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, None, None, 0)
# Z1_step2 = Step(Zephyr1, 0.45, 2, 1.05, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, material[6], fairing_shape[0], 0)
# Zephyr1_steps = [Z1_step1, Z1_step2]
# LVSteps.append(Zephyr1_steps)

# Zephyr2 = LaunchVehicle('Zephyr-2', 1.3, material[0], 3, [330, 380, 380], [0.11, 0.11, 0.11], 95, Two)
# LaunchVehicles.append(Zephyr2)
# #                    # (r, laststep, TW, engine,  dome shape,propellant, body material, tank material, number of engines, t_start, fairing_material, fairing_shape) 
# Z2_step1 = Step(Zephyr2, 0.45, 1, 1.4, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, None, None, 0)
# Z2_step2 = Step(Zephyr2, 0.45, 2, 1.05, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, None, None, 0)
# Z2_step3 = Step(Zephyr2, 0.39, 3, 0.9, engine[0], dome_shape[1], propellant[1], pressurant[0], ox_insul, [material[3], material[3]], 1, 1, 1, material[6], fairing_shape[0], 0)
# Zephyr2_steps = [Z2_step1, Z2_step2, Z2_step3]
# LVSteps.append(Zephyr2_steps)

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
#runTrajectory() # function runs trajectory for launch vehicles specified inside function
# for i in range(len(LaunchVehicles)):
#     LV = LaunchVehicles[i]
#     LV.massMoments(loads_conditions[1]) # wind loads


# Tkinter Tutorial
# Declare a window
root = Tk() # root window
# Create a Widget
myLabel1 = Label(root, text = "I'm on the first row")
myLabel2 = Label(root, text = "I'm on the second row")
# Grid positioning System
myLabel1.grid(row=0,column=0) # 'packs' the widget 'myLabel' into the window 'root' (shows widget on screen)
myLabel2.grid(row=0,column=1) # The column's width is the width of the widest widget in the row(=0 if there is no widget on previous columns)


# Create a Button
def myClick(): # Define the function to be executed after the button is pressed
    myLabeltext = Enter.get()
    myLabel = Label(root, text= myLabeltext)
    myLabel.grid(row=1,column=1)
myButton = Button(root, text="Enter your name", padx=50, pady=100, command=myClick, fg="white", bg = "black")
# 'padx' and 'pady' define button size 
# 'state=DISABLED' makes the button unclickable
# 'command=functionName' executes the function when user clicks on it.
# 'bg="color"' and 'fg="color"' changes the color of the background (button) and the foreground (text)
    # list of colors: red, blue, black, white, ..., and hex color codes: #000000, #ffffff
myButton.grid(row=1,column=0)

# Input Buttons (Entry Bars)
Enter = Entry(root, width=50, fg = "black", bg ="white", borderwidth = 5)
Enter.grid(row=2,column=0)
Enter.insert(0, "Enter Your Name:") # Insert default text inside entry bar

root.mainloop() # continuously iterates window
print('The program is finished!')