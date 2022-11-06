# LaunchVehicleDesign

The 'MainLVDesign.py' script is the main script used to generate mass budgets and size launch vehicles to complete a mission for a given altitude, inclination, and payload.

1. First the launch vehicle's propellant mass budgets are estimated by calling a MATLAB script in the function "generateLVMassEstimates(missions)", which takes a list of mission objects as the parameter.

2. Mission objects detail the orbital requirements and the required delta-v to complete the mission.

3. After the propellant mass budgets are calculated, the "LaunchVehicle" objects are initialized, which are the parents of the "Step" objects.

4. The "Step" objects initialize parameters such as propellant choice, gas used to pressurize tanks, materials, dome shape, and number of engines. The geometry of the step's body and tank components are then sized to accomodate the required volume of propellant.

5. After each "Step" is sized, the interstages of the "LaunchVehicle" are sized and the moments of inertia of each launch vehicle component are tabulated into a table.

6. A MATLAB script is called in "runTrajectory()" to run trajectory simulations of the "LaunchVehicle" objects using a gravity-turn and the optimal trajectory is selected based on quickest time to orbit and largest amount of propellant left for reserves.

7. The conditions at the point of maximum dynamic pressure (Max-Q) for the optimal trajectory is tabulated, such as: angle of attack, vehicle velocity, wind speed, left-over propellant mass.

# SCRIPT DOES NOT YET AUTOMATE:::
1. Excel is used to make the Axial Load, Shear Load, and Bending Moment diagrams from the output of the MatLab trajectory script. are made for the "LaunchVehicle" at the Max-Q condition, and Ground Wind Load Condition in Excel.  

2. Body component thickness selection: The stresses along the launch vehicle's height are manually tabulated in Excel, and the thicknesses of body components are found using the goal-seek function such that there is a safety margin of 1.5 for the ultimate compression strength of the material

3. The script also exports a PPT presentation slide of a to-scale side-view profile of the "LaunchVehicle". <-- Work still must be done to position PPT shapes correctly
