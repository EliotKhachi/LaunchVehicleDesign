%This Script is used to run the Trajectory Function and find and graph the
%optimal trajectory
clear all

%Vandenberg: Lat = 34.60 deg.
%Kennedy Space Center: Lat = 28.50 deg.

%Preliminary Input for Minerva 1: ('Minerva',[0.2; 1.0; 3282.7; 295.443; (1.4*9.80665*5165.5) ; 296.1] , [0.2; 1.0; 1852.8; 164.0526; (1.05*9.80665*1852.8); 359.1] , [0; 0; 0; 0; 0; 0], 1, 28.50)
%Preliminary Input for Minerva 2: ('Minerva',[0.2; 1.0; 3282.7; 295.443; (1.4*9.80665*5656.4) ; 296.1] , [0.2; 1.0; 1822.8; 164.0526; (1.05*9.80665*2373.6848); 359.1] , [0.2; 1.0; 550.8848; 41.0296; (0.9*9.80665*550.8848); 359.1] , 2, 34.60)
%Preliminary Input for Latona 1:  ('Latona', [0.2; 1.0; 1013.9; 60.832; (1.4*9.80665*1254.1) ; 265] , [0.2; 1.0; 240.2467; 18.9222;(1.05*9.80665*240.2467) ; 380] , [0; 0; 0; 0; 0; 0], 1, 28.50)
%Preliminary Input for Latona 2:  ('Latona', [0.2; 1.0; 3601.2; 216.0728; (1.4*9.80665*4920.3) ; 265] , [0.2; 1.0; 1013.9; 60.8320; (1.05*9.80665*1463.0879) ; 265] , [0.2; 1.0; 305.2497; 18.9222;(0.9*9.80665*305.2497) ; 380] , 2, 34.60)
% Please Note: It is assumed that for Latona 2 there are
% 4 SRBs

%Output as X by 10 double that is columns for the following:
%scale_factor_1, scale_factor_2, scale_factor_3, pitch_kick, mleft_1, mleft_2, mleft_3, count, delta_v_total, delta_v_circularization

%-------- Variables for Inputs for Function defined here (in kg,m,s format)
%% 

%% 
One_dv_params = readtable('LVMasses\OnedVParameters.csv'); % load Mission 1 dv requirements csv into table
Two_dv_params = readtable('LVMasses\TwodVParameters.csv'); % load Mission 2 dv requirements csv into table 
Latona1_params = readtable('LVTrajectory\Latona-1TrajReqs.csv');
Latona2_params = readtable('LVTrajectory\Latona-2TrajReqs.csv');
Minerva1_params = readtable('LVTrajectory\Minerva-1TrajReqs.csv');
Minerva2_params = readtable('LVTrajectory\Minerva-2TrajReqs.csv');
%--Minerva Inputs--
grav_sea_level = 9.80665;
Mass_Payload_1 = One_dv_params{1,9};
Mass_Payload_2 = Two_dv_params{1,9};

Minerva_1_Cd2 = 0.2;
Minerva_1_Radius2 = Minerva1_params{2,2}; 
Minerva_1_mi2 = Minerva1_params{2,3};
Minerva_1_mf2 = Minerva1_params{2,4} + Mass_Payload_1;
Minerva_1_mprop2 = Minerva_1_mi2 - Minerva_1_mf2;
Minerva_1_T2 = 1.05 * (Minerva_1_mi2) * grav_sea_level;
Minerva_1_Isp2 = Minerva1_params{2,8};

Minerva_1_Cd1 = 0.2;
Minerva_1_Radius1 = Minerva1_params{1,2};
Minerva_1_mi1 = Minerva1_params{1,3};
Minerva_1_mf1 = Minerva1_params{1,4};
Minerva_1_mprop1 = Minerva_1_mi1 - Minerva_1_mf1;
Minerva_1_T1 = 1.4 * (Minerva_1_mi1 + Minerva_1_mi2) * grav_sea_level;
Minerva_1_Isp1 = Minerva1_params{1,8};

%---

Minerva_2_Cd3 = 0.2;
Minerva_2_Radius3 = Minerva2_params{3,2};
Minerva_2_mi3 = Minerva2_params{3,3};
Minerva_2_mf3 = Minerva2_params{3,4} + Mass_Payload_2;
Minerva_2_mprop3 = Minerva_2_mi3 - Minerva_2_mf3;
Minerva_2_T3 = 0.9 * (Minerva_2_mi3) * grav_sea_level;
Minerva_2_Isp3 = Minerva2_params{3,8};

Minerva_2_Cd2 = 0.2;
Minerva_2_Radius2 = Minerva2_params{2,2};
Minerva_2_mi2 = Minerva2_params{2,3};
Minerva_2_mf2 = Minerva2_params{2,4};
Minerva_2_mprop2 = Minerva_2_mi2 - Minerva_2_mf2;
Minerva_2_T2 = 1.05 * (Minerva_2_mi2 + Minerva_2_mi3) * grav_sea_level;
Minerva_2_Isp2 = Minerva2_params{2,8};

Minerva_2_Cd1 = 0.2;
Minerva_2_Radius1 = Minerva2_params{1,2};
Minerva_2_mi1 = Minerva2_params{1,3};
Minerva_2_mf1 = Minerva2_params{1,4};
Minerva_2_mprop1 = Minerva_2_mi1 - Minerva_2_mf1;
Minerva_2_T1 = 1.4 * (Minerva_2_mi1 + Minerva_2_mi2 + Minerva_2_mi3) * grav_sea_level;
Minerva_2_Isp1 = Minerva2_params{1,8};

%--Latona Inputs--

Latona_1_Cd2 = 0.2;
Latona_1_Radius2 = Latona1_params{2,2};
Latona_1_mi2 = Latona1_params{2,3};
Latona_1_mf2 = Latona1_params{2,4} + Mass_Payload_1;
Latona_1_mprop2 = Latona_1_mi2 - Latona_1_mf2;
Latona_1_T2 = 1.05 * (Latona_1_mi2) * grav_sea_level;
Latona_1_Isp2 = Latona1_params{2,8};

Latona_1_Cd1 = 0.2;
Latona_1_Radius1 = Latona1_params{1,2};
Latona_1_mi1 = Latona1_params{1,3};
Latona_1_mf1 = Latona1_params{1,4};
Latona_1_mprop1 = Latona_1_mi1 - Latona_1_mf1;
Latona_1_T1 = 1.4 * (Latona_1_mi1 + Latona_1_mi2) * grav_sea_level;
Latona_1_Isp1 = Latona1_params{1,8};

%--- %Latona 2 Step 1 is actually step '0'

Latona_2_Cd3 = 0.2;                   
Latona_2_Radius3 = Latona2_params{3,2};
Latona_2_mi3 = Latona2_params{3,3};
Latona_2_mf3 = Latona2_params{3,4} + Mass_Payload_2;
Latona_2_mprop3 = Latona_2_mi3 - Latona_2_mf3;
Latona_2_T3 = 0.9 * (Latona_2_mi3) * grav_sea_level;
Latona_2_Isp3 = Latona2_params{3,8};

Latona_2_Cd2 = 0.2;                   
Latona_2_Radius2 = Latona2_params{2,2};
Latona_2_mi2 = Latona2_params{2,3};
Latona_2_mf2 = Latona2_params{2,4};
Latona_2_mprop2 = Latona_2_mi2 - Latona_2_mf2;
Latona_2_T2 = Latona_1_T1;                      %Same solid rocket is used for mission 2 as for mission 1
Latona_2_Isp2 = Latona2_params{2,8};

Latona_2_Cd1 = 0.2;                   
Latona_2_Radius1 = Latona2_params{1,2};
Latona_2_mi1 = Latona2_params{1,3};
Latona_2_mf1 = Latona2_params{1,4};
Latona_2_mprop1 = Latona_2_mi1 - Latona_2_mf1;
Latona_2_T1 = 1.4 * (Latona_2_mi3 + Latona_2_mi2 + Latona_2_mi1) * grav_sea_level;
Latona_2_Main_Thrust = Latona_1_T1;
Latona_2_Booster_Thrust = Latona_2_T1 - Latona_2_Main_Thrust;
Latona_2_Isp1 = Latona2_params{1,8};
%----Using Trajectory Function with inputvariables above

Minerva_1_Trajectories = Trajectory_TL_0320('Minerva',[Minerva_1_Cd1; Minerva_1_Radius1; Minerva_1_mi1; Minerva_1_mf1; Minerva_1_T1 ; Minerva_1_Isp1] , [Minerva_1_Cd2; Minerva_1_Radius2; Minerva_1_mi2; Minerva_1_mf2; Minerva_1_T2; Minerva_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, 28.50);
Minerva_2_Trajectories = Trajectory_TL_0320('Minerva',[Minerva_2_Cd1; Minerva_2_Radius1; Minerva_2_mi1; Minerva_2_mf1; Minerva_2_T1 ; Minerva_2_Isp1] , [Minerva_2_Cd2; Minerva_2_Radius2; Minerva_2_mi2; Minerva_2_mf2; Minerva_2_T2; Minerva_2_Isp2] , [Minerva_2_Cd3; Minerva_2_Radius3; Minerva_2_mi3; Minerva_2_mf3; Minerva_2_T3; Minerva_2_Isp3] , 2, 34.60);
Latona_1_Trajectories = Trajectory_TL_0320('Latona', [Latona_1_Cd1; Latona_1_Radius1; Latona_1_mi1; Latona_1_mf1; Latona_1_T1 ; Latona_1_Isp1] , [Latona_1_Cd2; Latona_1_Radius2; Latona_1_mi2; Latona_1_mf2; Latona_1_T2; Latona_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, 28.50);
Latona_2_Trajectories = Trajectory_TL_0320('Latona', [Latona_2_Cd1; Latona_2_Radius1; Latona_2_mi1; Latona_2_mf1; Latona_2_T1 ; Latona_2_Isp1] , [Latona_2_Cd2; Latona_2_Radius2; Latona_2_mi2; Latona_2_mf2; Latona_2_T2; Latona_2_Isp2] , [Latona_2_Cd3; Latona_2_Radius3; Latona_2_mi3; Latona_2_mf3; Latona_2_T3; Latona_2_Isp3] , 2, 34.60);

%Sort_time_Minerva_1 = sortrows(Minerva_1_Trajectories, [8 10]); 
%Sort_time_Minerva_2 = sortrows(Minerva_2_Trajectories, [8 10]);
%Sort_time_Latona_1 = sortrows(Latona_1_Trajectories, [8 10]);
%Sort_time_Latona_2 = sortrows(Latona_2_Trajectories, [8 10]);

%Sort_DVcirc_Minerva_1 = sortrows(Minerva_1_Trajectories, [10 8]); 
%Sort_DVcirc_Minerva_2 = sortrows(Minerva_2_Trajectories, [10 8]);
%Sort_DVcirc_Latona_1 = sortrows(Latona_1_Trajectories, [10 8]);
%Sort_DVcirc_Latona_2 = sortrows(Latona_2_Trajectories, [10 8]);

[Minimum_Values_1 , Index_1] = min(Minerva_1_Trajectories, [], 1);
[Minimum_Values_2 , Index_2] = min(Minerva_2_Trajectories, [], 1);
[Minimum_Values_3 , Index_3] = min(Latona_1_Trajectories, [], 1);
[Minimum_Values_4 , Index_4] = min(Latona_2_Trajectories, [], 1);

Optimal_Minerva_1_Result = Minerva_1_Trajectories(Index_1(10),:)
Optimal_Minerva_2_Result = Minerva_2_Trajectories(Index_2(10),:)
Optimal_Latona_1_Result = Latona_1_Trajectories(Index_3(10),:)
Optimal_Latona_2_Result = Latona_2_Trajectories(Index_4(10),:)

%-------------------------------------------------------------------------
%This code runs the Trajectory function once for when the optimal
%trajectory is found and graphs it

Optimal_Minerva_1_Trajectory = Optimal_Trajectory('Minerva',[Minerva_1_Cd1; Minerva_1_Radius1; Minerva_1_mi1; Minerva_1_mf1; Minerva_1_T1 ; Minerva_1_Isp1] , [Minerva_1_Cd2; Minerva_1_Radius2; Minerva_1_mi2; Minerva_1_mf2; Minerva_1_T2; Minerva_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, 28.50,Optimal_Minerva_1_Result);
Optimal_Minerva_2_Trajectory = Optimal_Trajectory('Minerva',[Minerva_2_Cd1; Minerva_2_Radius1; Minerva_2_mi1; Minerva_2_mf1; Minerva_2_T1 ; Minerva_2_Isp1] , [Minerva_2_Cd2; Minerva_2_Radius2; Minerva_2_mi2; Minerva_2_mf2; Minerva_2_T2; Minerva_2_Isp2] , [Minerva_2_Cd3; Minerva_2_Radius3; Minerva_2_mi3; Minerva_2_mf3; Minerva_2_T3; Minerva_2_Isp3] , 2, 34.60,Optimal_Minerva_2_Result);
Optimal_Latona_1_Trajectory = Optimal_Trajectory('Latona', [Latona_1_Cd1; Latona_1_Radius1; Latona_1_mi1; Latona_1_mf1; Latona_1_T1 ; Latona_1_Isp1] , [Latona_1_Cd2; Latona_1_Radius2; Latona_1_mi2; Latona_1_mf2; Latona_1_T2; Latona_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, 28.50, Optimal_Latona_1_Result);
Optimal_Latona_2_Trajectory = Optimal_Trajectory('Latona', [Latona_2_Cd1; Latona_2_Radius1; Latona_2_mi1; Latona_2_mf1; Latona_2_T1 ; Latona_2_Isp1] , [Latona_2_Cd2; Latona_2_Radius2; Latona_2_mi2; Latona_2_mf2; Latona_2_T2; Latona_2_Isp2] , [Latona_2_Cd3; Latona_2_Radius3; Latona_2_mi3; Latona_2_mf3; Latona_2_T3; Latona_2_Isp3] , 2, 34.60,Optimal_Latona_2_Result);

Optimal_Minerva_1_Trajectory = Optimal_Minerva_1_Trajectory(1:Minimum_Values_1(8), :);
Optimal_Minerva_2_Trajectory = Optimal_Minerva_2_Trajectory(1:Minimum_Values_2(8), :);
Optimal_Latona_1_Trajectory = Optimal_Latona_1_Trajectory(1:Minimum_Values_3(8), :);
Optimal_Latona_2_Trajectory = Optimal_Latona_2_Trajectory(1:Minimum_Values_4(8), :);


scales1 = Optimal_Minerva_1_Result(1:7);
scales2 = Optimal_Minerva_2_Result(1:7);
scales3 = Optimal_Latona_1_Result(1:7);
scales4 = Optimal_Latona_2_Result(1:7);

plotTraj('Minerva',[Minerva_1_Cd1; Minerva_1_Radius1; Minerva_1_mi1; Minerva_1_mf1; Minerva_1_T1 ; Minerva_1_Isp1] , [Minerva_1_Cd2; Minerva_1_Radius2; Minerva_1_mi2; Minerva_1_mf2; Minerva_1_T2; Minerva_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, scales1);
plotTraj('Minerva',[Minerva_2_Cd1; Minerva_2_Radius1; Minerva_2_mi1; Minerva_2_mf1; Minerva_2_T1 ; Minerva_2_Isp1] , [Minerva_2_Cd2; Minerva_2_Radius2; Minerva_2_mi2; Minerva_2_mf2; Minerva_2_T2; Minerva_2_Isp2] , [Minerva_2_Cd3; Minerva_2_Radius3; Minerva_2_mi3; Minerva_2_mf3; Minerva_2_T3; Minerva_2_Isp3] , 2, scales2);
plotTraj('Latona', [Latona_1_Cd1; Latona_1_Radius1; Latona_1_mi1; Latona_1_mf1; Latona_1_T1 ; Latona_1_Isp1] , [Latona_1_Cd2; Latona_1_Radius2; Latona_1_mi2; Latona_1_mf2; Latona_1_T2; Latona_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, scales3);
plotTraj('Latona', [Latona_2_Cd1; Latona_2_Radius1; Latona_2_mi1; Latona_2_mf1; Latona_2_T1 ; Latona_2_Isp1] , [Latona_2_Cd2; Latona_2_Radius2; Latona_2_mi2; Latona_2_mf2; Latona_2_T2; Latona_2_Isp2] , [Latona_2_Cd3; Latona_2_Radius3; Latona_2_mi3; Latona_2_mf3; Latona_2_T3; Latona_2_Isp3] , 2, scales4);




% ----------------------- Plotting --------------------------------------

% t_Minerva_1 = Optimal_Minerva_1_Trajectory(:,1);
% v_Minerva_1 = Optimal_Minerva_1_Trajectory(:,2);
% m_Minerva_1 = Optimal_Minerva_1_Trajectory(:,3);
% h_Minerva_1 = Optimal_Minerva_1_Trajectory(:,4);
% x_Minerva_1 = Optimal_Minerva_1_Trajectory(:,5);
% gamma_Minerva_1 = Optimal_Minerva_1_Trajectory(:,6);
% a_Minerva_1 = Optimal_Minerva_1_Trajectory(:,7);
% q_Minerva_1 = Optimal_Minerva_1_Trajectory(:,8);

% yyaxis left
% plot(t_Minerva_1, rad2deg(gamma_Minerva_1), t_Minerva_1, h_Minerva_1/1000, t_Minerva_1, a_Minerva_1)
% xlabel('time (s)')
% ylabel('Altitude (km), Flight Path Angle (Deg), Acceleration (m/s^2)')
% ylim([0 inf])
% 
% yyaxis right
% plot(t_Minerva_1, v_Minerva_1, t_Minerva_1, q_Minerva_1/1000)
% ylabel('Velocity (m/s), dynamic pressure (kPa)')
% legend('gamma', 'altitude', 'acceleration', 'velocity', 'dynamic pressure')
% %title([name ' ' string(mission) ' Trajectory'])
% ylim([0 inf])

% t_Latona_1 = Optimal_Latona_1_Trajectory(:,1);
% v_Latona_1 = Optimal_Latona_1_Trajectory(:,2);
% m_Latona_1 = Optimal_Latona_1_Trajectory(:,3);
% h_Latona_1 = Optimal_Latona_1_Trajectory(:,4);
% x_Latona_1 = Optimal_Latona_1_Trajectory(:,5);
% gamma_Latona_1 = Optimal_Latona_1_Trajectory(:,6);
% a_Latona_1 = Optimal_Latona_1_Trajectory(:,7);
% q_Latona_1 = Optimal_Latona_1_Trajectory(:,8);

% ------

% yyaxis left
% plot(t_Latona_1, rad2deg(gamma_Latona_1), t_Latona_1, h_Latona_1/1000, t_Latona_1, a_Latona_1)
% xlabel('time (s)')
% ylabel('Altitude (km), Flight Path Angle (Deg), Acceleration (m/s^2)')
% ylim([0 inf])
% 
% yyaxis right
% plot(t_Latona_1, v_Latona_1, t_Latona_1, q_Latona_1/1000)
% ylabel('Velocity (m/s), dynamic pressure (kPa)')
% legend('gamma', 'altitude', 'acceleration', 'velocity', 'dynamic pressure')
% %title([name ' ' string(mission) ' Trajectory'])
% ylim([0 inf])
% grid on
% xline(0)

%-----------Imbeded Function ---------------------------------------------
    function Optimal_Result = Optimal_Trajectory(name, step1, step2, step3, mission, launch_latitude, Optimals)
%% Earth and Launch Site Inputs
    %% Earth and Launch Site Inputs
    mu = 3.986e14; %m^3/s^2             Gravitational Parameter of Earth
    g0 = 9.80665; %m/s^2                Local Gravity at S.L.
    R_earth = 6378000; %m               Radius of Earth
    h0 = 7640; %m                       Scale Height (Wikipedia: https://en.wikipedia.org/wiki/Scale_height)
    L0 = launch_latitude; %deg, N       Latitude of Launch Site
    rho0 = 1.225; %kg/m^3               Air Density at S.L.
    v_ls = 465.1*cos(deg2rad(L0)); %m/s Speed of Launch Site
    
    if mission == 1
        final_alt = 500000; %m          Final Orbit
        rf = final_alt+R_earth; %m      Distance btwn Final Orbit Alt and Center of Earth
        inc = 60; %deg                  Inclination
    elseif mission == 2
        final_alt = 550000; %m          Final Orbit
        rf = final_alt+R_earth; %m      Distance btwn Final Orbit Alt and Center of Earth
        inc = 95; %deg                  Inclination
    else
        warning('This mission does not exist in this simulation!')
    end
    
    %% Iterators
    stage1_Thrust_Scaling_Factor = .9:.1:1.1;
    %stage1_Thrust_Scaling_Factor = 0.9;
    
    stage2_Thrust_Scaling_Factor = .6:.15:1.1;
    %stage2_Thrust_Scaling_Factor = 0.6;
    
    stage3_Thrust_Scaling_Factor = .5:.15:1.1;
    %stage3_Thrust_Scaling_Factor = 0.5;
    
    pitchkick_spread = 0.01:.005:0.1;
    pitchkick_spread = deg2rad(pitchkick_spread);
        

    stage1_mleft = 0;   %percentage                            %Stage 1 Percentage of fuel remaining
    
    
    if mission == 2
        stage2_mleft = 0; %percentage                          %Stage 2 Percentage of fuel remaining
        stage3_mleft = 0.1:.05:.4; %percentage                 %Stage 3 Percentage of fuel remaining                   
    else
        stage2_mleft = 0.1:.05:.4; %percentage                 %Stage 2 Percentage of fuel remaining
        stage3_mleft = 0; %percentage                          %Stage 3 Percentage of fuel remaining                   
    end

    %% Initial Conditions
    t0 = 0; %s                                          %Initial Time
    dt = 1; %s                                          %Change in Time    
    gamma0 = pi/2; %rad                                 %Initial Flight Path Angle
    gamma_dot0 = 0; %rad/s                              %Initial Change in Flight Path Angle
    drag0 = 0; %N                                       %Initial Drag
    x0 = 0; %m                                          %Initial Downrange Distance
    xdot0 = 0; %m/s                                     %Initial Downrange Speed
    h00 = 0; %m                                         %Initial Vertical Distance (Launch Site Altitude)
    hdot0 = 0; %m/s                                     %Initial Vertical Speed 
    v0 = 0; %m/s                                        %Initial Speed
    q0 = 0; %Pa                                        %Initial Dynamic Pressure
    %Label_Results = ['Stage 1 Throttle', 'Stage 2 Throttle' , 'Stage 3 Throttle', 'Pitch-Kick' , 'Check', 'Delta-V to Circularize'];
    %Results = (Label_Results);
    Results = [];
    Record = [];
    %% Simulation
    
     scale_factor_1 = Optimals(1);
        scale_factor_2 = Optimals(2);
             scale_factor_3 = Optimals(3);
                 pitch_kick = Optimals(4);
                     mleft_1 = Optimals(5);
                         mleft_2 = Optimals(6);
                             mleft_3 = Optimals(7);
                                %% Array Initialization                              
                                t = zeros(600, 1);
                                t(1) = t0;
                                gamma_dot = zeros(600, 1);
                                gamma_dot(1) = gamma_dot0;
                                gamma = zeros(600, 1);
                                gamma(1) = gamma0;
                                xdot = zeros(600, 1);
                                xdot(1) = xdot0;
                                x = zeros(600, 1);
                                x(1) = x0;
                                hdot = zeros(600, 1);
                                hdot(1) = hdot0;
                                h = zeros(600, 1);
                                h(1) = h00;
                                rho = zeros(600, 1);
                                rho(1) = rho0;
                                drag = zeros(600, 1);
                                drag(1) = drag0;
                                g = zeros(600, 1);
                                g(1) = g0;
                                v = zeros(600, 1);
                                v(1) = v0;
                                q = zeros(600, 1);
                                q(1) = q0;
                                
                                %% LV Configuration Data
                                stage1_Cd = step1(1);                                    %Stage 1 Coefficient of Drag
                                stage1_Radius = step1(2); %m                             %Stage 1 Radius
                                stage1_mi = step1(3)+step2(3)+step3(3); %kg              %Stage 1 Initial Mass (Stage l Struct, Stage 1 Fuel, Stage 2 Struct, Stage 2 Fuel)
                                if mission == 2 && strcmp(name, 'Latona')
                                   stage1_S = pi*(stage1_Radius^2)*4 + pi*step2(2)^2; %m^2 %Stage 1 Cross-Sectional Area
                                else
                                   stage1_S = pi*stage1_Radius^2; %m^2                   %Stage 1 Cross-Sectional Area
                                end
                                stage1_Thrust = scale_factor_1*step1(5); %N              %Stage 1 Thrust
                                stage1_Isp = step1(6); %s                                %Stage 1 Isp
                                
                                if mission == 2 && strcmp(name, 'Latona')
                                    stage1_Main_Thrust = 1.4 * (step2(3) + step3(3));
                                    stage1_Main_mdot = stage1_Main_Thrust/(step2(6)*g0);
                                    stage1_Booster_Thrust_tot = step1(5) - stage1_Main_Thrust;
                                    stage1_mdot = stage1_Booster_Thrust_tot/(stage1_Isp*g0);
                                    stage1_Booster_mdot_indv = stage1_mdot/4;
                                else
                                    stage1_mdot = stage1_Thrust/(stage1_Isp*g0); %kg/s   %Stage 1 Mass Flow
                                end
                                if mission == 2 && strcmp(name, 'Latona')
                                    s1prop = step1(3) - step1(4);
                                    s1burntime = s1prop/stage1_mdot;
                                    s2prop = step2(3) - step2(4);
                                    s2prop_burn = stage1_Main_mdot * s1burntime;
                                    stage1_mf = step1(4) + (step2(3) - s2prop_burn)+step3(3);
                                else
                                    stage1_mf = step1(4)+step2(3)+step3(3); %kg          %Stage 1 Final Mass (Stage 1 Struct, Stage 2 Fuel, Stage 2 Struct)     
                                end
                                
                                stage2_Cd = step2(1);                                    %Stage 2 Coefficient of Drag
                                stage2_Radius = step2(2); %m                             %Stage 2 Radius
                                if mission == 2 && strcmp(name, 'Latona')
                                    stage2_mi = (step2(3) - s2prop_burn)+step3(3);
                                else
                                    stage2_mi = step2(3)+step3(3); %kg                   %Stage 2 Initial Mass (Stage l Struct, Stage 1 Fuel, Stage 2 Struct, Stage 2 Fuel)
                                end
                                stage2_mf = step2(4)+step3(3); %kg                       %Stage 2 Final Mass (Stage 1 Struct, Stage 2 Fuel, Stage 2 Struct)
                                stage2_S = pi*stage2_Radius^2; %m^2                      %Stage 2 Cross-Sectional Area
                                stage2_Thrust = scale_factor_2*step2(5); %N                          %Stage 2 Thrust
                                stage2_Isp = step2(6); %s                                %Stage 2 Isp
                                if mission == 2 && strcmp(name, 'Latona')
                                    stage2_mdot = stage1_Main_mdot;
                                else
                                    stage2_mdot = stage2_Thrust/(stage2_Isp*g0); %kg/s   %Stage 2 Mass Flow
                                end
                                
                                stage3_Cd = step3(1);                                    %Stage 3 Coefficient of Drag
                                stage3_Radius = step3(2); %m                             %Stage 3 Radius
                                stage3_mi = step3(3); %kg                                %Stage 3 Initial Mass (Stage l Struct, Stage 1 Fuel, Stage 2 Struct, Stage 2 Fuel)
                                stage3_mf = step3(4); %kg                                %Stage 3 Final Mass (Stage 1 Struct, Stage 2 Fuel, Stage 2 Struct)
                                stage3_S = pi*stage3_Radius^2; %m^2                      %Stage 3 Cross-Sectional Area
                                stage3_Thrust = scale_factor_3*step3(5); %N                           %Stage 3 Thrust
                                stage3_Isp = step3(6); %s                                %Stage 3 Isp
                                if stage3_mi ~= 0
                                    stage3_mdot = stage3_Thrust/(stage3_Isp*g0); %kg/s   %Stage 3 Mass Flow
                                else
                                    stage3_mdot = 0; %kg/s                               %Stage 3 Mass Flow
                                end
                                
                                a0 = (stage1_Thrust/stage1_mi) - g0*sin(gamma0); %m/s^2  %Initial Acceleration
                                a = zeros(600, 1);
                                a(1) = a0;
                                m = zeros(600, 1);
                                m(1) = stage1_mi;
                                
                                count = 1;                                               %Begin Count Here
                                while gamma(count) > deg2rad(1)
                                    t(count+1) = t(count)+1;
                                    if m(count) > stage1_mf+(mleft_1*(stage1_mi-stage1_mf))%The Order of calculating these terms are important
                                        v(count+1) = v(count) + a(count)*dt;
                                        if mission == 2 && strcmp(name, 'Latona')
                                           m(count+1) = m(count) - stage1_mdot*dt - stage2_mdot*dt;
                                        else
                                           m(count+1) = m(count) - stage1_mdot*dt; 
                                        end
                                        h(count+1) = h(count) + hdot(count)*dt;
                                        x(count+1) = x(count) + xdot(count)*dt;
                                        if h(count) >= 400 && h(count) <= 600
                                            gamma(count+1) = -pitch_kick + gamma(count) + gamma_dot(count)*dt;
                                        else
                                            gamma(count+1) = gamma(count) + gamma_dot(count)*dt;
                                        end
                                        rho(count+1) = rho0*exp(-h(count+1)/h0);
                                        drag(count+1) = .5*stage1_Cd*stage1_S*rho(count+1)*v(count+1)^2;
                                        g(count+1) = g0/((1+(h(count+1)/R_earth))^2);
                                        gamma_dot(count+1) = -((g(count+1)/v(count+1)) - (v(count+1)/(R_earth+h(count+1))))*cos(gamma(count+1));
                                        a(count+1) = stage1_Thrust/m(count+1) - drag(count+1)/m(count+1) - (g(count+1)*sin(gamma(count+1)));
                                        xdot(count+1) = v(count+1)*cos(gamma(count+1))*R_earth/(R_earth+h(count+1));
                                        hdot(count+1) = v(count+1)*sin(gamma(count+1));
                                        q(count+1) = drag(count+1)/(stage1_Cd*stage1_S);

                                    elseif m(count) > stage2_mi && m(count) <= stage1_mf+(mleft_1*(stage1_mi-stage1_mf))
                                        ttemp = t(count);    
                                        for j = 1:5
                                            v(count+1) = v(count) + a(count)*dt;
                                            m(count+1) = stage2_mi;
                                            h(count+1) = h(count) + hdot(count)*dt;
                                            x(count+1) = x(count) + xdot(count)*dt;
                                            gamma(count+1) = gamma(count) + gamma_dot(count)*dt;
                                            rho(count+1) = rho0*exp(-h(count+1)/h0);
                                            drag(count+1) = .5*stage2_Cd*stage2_S*rho(count+1)*v(count+1)^2;
                                            g(count+1) = g0/((1+(h(count+1)/R_earth))^2);
                                            gamma_dot(count+1) = -((g(count+1)/v(count+1))-(v(count+1)/(R_earth+h(count+1))))*cos(gamma(count+1));
                                            a(count+1) = (-drag(count+1))/m(count+1) - (g(count+1)*sin(gamma(count+1)));
                                            xdot(count+1) = v(count+1)*cos(gamma(count+1))*R_earth/(R_earth+h(count+1));
                                            hdot(count+1) = v(count+1)*sin(gamma(count+1));
                                            q(count+1) = drag(count)/(stage2_Cd*stage2_S);
                                            t(count+1) = ttemp + j;
                                        end

                                    elseif m(count) > stage2_mf+(mleft_2*(stage2_mi-stage2_mf))
                                        v(count+1) = v(count) + a(count)*dt;
                                        m(count+1) = m(count) - stage2_mdot*dt;
                                        h(count+1) = h(count) + hdot(count)*dt;
                                        x(count+1) = x(count) + xdot(count)*dt;
                                        gamma(count+1) = gamma(count) + gamma_dot(count)*dt;
                                        rho(count+1) = rho0*exp(-h(count+1)/h0);
                                        drag(count+1) = .5*stage2_Cd*stage2_S*rho(count+1)*v(count+1)^2;
                                        g(count+1) = g0/((1+(h(count+1)/R_earth))^2);
                                        gamma_dot(count+1) = -((g(count+1)/v(count+1))-(v(count+1)/(R_earth+h(count+1))))*cos(gamma(count+1));
                                        a(count+1) = (stage2_Thrust-drag(count+1))/m(count+1) - (g(count+1)*sin(gamma(count+1)));
                                        xdot(count+1) = v(count+1)*cos(gamma(count+1))*R_earth/(R_earth+h(count+1));
                                        hdot(count+1) = v(count+1)*sin(gamma(count+1)); 
                                        q(count+1) = drag(count+1)/(stage2_Cd*stage2_S);

                                    elseif m(count) > stage3_mi && m(count) <= stage2_mf && stage3_mi ~= 0
                                        ttemp = t(count);    
                                        for j = 1:5
                                            v(count+1) = v(count) + a(count)*dt;
                                            m(count+1) = stage3_mi;
                                            h(count+1) = h(count) + hdot(count)*dt;
                                            x(count+1) = x(count) + xdot(count)*dt;
                                            gamma(count+1) = gamma(count) + gamma_dot(count)*dt;
                                            rho(count+1) = rho0*exp(-h(count+1)/h0);
                                            drag(count+1) = .5*stage3_Cd*stage3_S*rho(count+1)*v(count+1)^2;
                                            g(count+1) = g0/((1+(h(count+1)/R_earth))^2);
                                            gamma_dot(count+1) = -((g(count+1)/v(count+1))-(v(count+1)/(R_earth+h(count+1))))*cos(gamma(count+1));
                                            a(count+1) = (-drag(count+1))/m(count+1) - (g(count+1)*sin(gamma(count+1)));
                                            xdot(count+1) = v(count+1)*cos(gamma(count+1))*R_earth/(R_earth+h(count+1));
                                            hdot(count+1) = v(count+1)*sin(gamma(count+1));
                                            q(count+1) = drag(count+1)/(stage3_Cd*stage3_S);
                                            t(count+1) = ttemp + j;
                                        end

                                    elseif m(count) > stage3_mf+(mleft_3*(stage3_mi-stage3_mf)) && stage3_mi ~= 0
                                        v(count+1) = v(count) + a(count)*dt;
                                        m(count+1) = m(count) - stage3_mdot*dt;
                                        h(count+1) = h(count) + hdot(count)*dt;
                                        x(count+1) = x(count) + xdot(count)*dt;
                                        gamma(count+1) = gamma(count) + gamma_dot(count)*dt;
                                        rho(count+1) = rho0*exp(-h(count+1)/h0);
                                        drag(count+1) = .5*stage3_Cd*stage3_S*rho(count+1)*v(count+1)^2;
                                        g(count+1) = g0/((1+(h(count+1)/R_earth))^2);
                                        gamma_dot(count+1) = -((g(count+1)/v(count+1))-(v(count+1)/(R_earth+h(count+1))))*cos(gamma(count+1));
                                        a(count+1) = (stage3_Thrust-drag(count+1))/m(count+1) - (g(count+1)*sin(gamma(count+1)));
                                        xdot(count+1) = v(count+1)*cos(gamma(count+1))*R_earth/(R_earth+h(count+1));
                                        hdot(count+1) = v(count+1)*sin(gamma(count+1));
                                        q(count+1) = drag(count)/(stage3_Cd*stage3_S);

                                    else
                                        v(count+1) = v(count) + a(count)*dt;
                                        m(count+1) = m(count);
                                        h(count+1) = h(count) + hdot(count)*dt;
                                        x(count+1) = x(count) + xdot(count)*dt;
                                        gamma(count+1) = gamma(count) + gamma_dot(count)*dt;
                                        rho(count+1) = rho0*exp(-h(count+1)/h0);
                                        drag(count+1) = 0;
                                        g(count+1) = g0/((1+(h(count+1)/R_earth))^2);
                                        gamma_dot(count+1) = -((g(count+1)/v(count+1))-(v(count+1)/(R_earth+h(count+1))))*cos(gamma(count+1));
                                        a(count+1) = -(g(count+1)*sin(gamma(count+1)));
                                        xdot(count+1) = v(count+1)*cos(gamma(count+1))*R_earth/(R_earth+h(count+1));
                                        hdot(count+1) = v(count+1)*sin(gamma(count+1));
                                        q(count+1) = 0;
                                    end
                                    if count >= 3000
                                    break
                                    else
                                        count = count + 1;
                                    end
                                end
%                                 if length(v) ~= length(t)
%                                     t(count) = [];
%                                     if length(v) ~= length(t)
%                                         t(count) = [];
%                                         if length(v) ~= length(t)
%                                             warning('Time and Vector Length disagree. Please fix')
%                                         end
%                                     end
%                                 end


                                
                                %% Hohmann Transfer
                                r_periapsis = h(count) + R_earth;
                                r_apoapsis = rf;
                                eccentricity_transfer = (r_apoapsis - r_periapsis)/(r_apoapsis + r_periapsis);

%                                 if v(count) > 5000
%                                     disp('We are going too fast at this time. Please try again.')
%                                 end
                                
                                if inc < 90
                                    delta_v_circularization = sqrt(mu/r_periapsis) - (v(count) + v_ls);
                                else
                                    delta_v_circularization = sqrt(mu/r_periapsis) - (v(count) - v_ls);
                                end
                                delta_v_1 = sqrt(mu/r_periapsis) * (sqrt((2*r_apoapsis)/(r_periapsis+r_apoapsis))-1);
                                delta_v_2 = sqrt(mu/r_apoapsis) * (1-sqrt((2*r_periapsis)/(r_periapsis+r_apoapsis)));

                                delta_v_total = delta_v_1 + delta_v_2 + delta_v_circularization;

%                                 if stage3_mi == 0
%                                     MR_H_Transfer = exp(delta_v_total / (g0*stage2_Isp));
%                                     if MR_H_Transfer > ((stage2_mf+mleft_2*(stage2_mi-stage2_mf))/(stage2_mf))
%                                         check = 0;
%                                     else
%                                         check = 1;
%                                     end
%                                 else
%                                     MR_H_Transfer = exp(delta_v_total / (g0*stage3_Isp));
%                                     if MR_H_Transfer > ((stage3_mf+mleft_3*(stage3_mi-stage3_mf))/(stage3_mf))
%                                         check = 0;
%                                     else
%                                         check = 1;
%                                     end
%                                 end
%                                 
%                                 if h(count) >= 250000
%                                    check = 0; 
%                                 end
%                                 
%                                 if delta_v_total <= 0
%                                     check = 0;
%                                 end
%                                 
%                                 if check == 1
%                                     Results = [Results; scale_factor_1, scale_factor_2, scale_factor_3, pitch_kick, mleft_1, mleft_2, mleft_3, count, delta_v_total, delta_v_circularization];
%                                 end
                                %Record = [Record; scale_factor_1, scale_factor_2, scale_factor_3, pitch_kick, mleft_1, mleft_2, mleft_3, check, count, delta_v_total];
                                Optimal_Result = [t,v,m,h,x,gamma,a,q];
                                %clear a v gamma_dot gamma xdot x hdot h rho drag g m q
    end

