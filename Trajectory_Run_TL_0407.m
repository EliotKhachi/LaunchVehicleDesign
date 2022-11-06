%This Script is used to run the Trajectory Function and find and graph the
%optimal trajectory
clear all

%Vandenberg: Lat = 34.60 deg.
%Kennedy Space Center: Lat = 28.50 deg.

%Output as X by 10 double that is columns for the following:
%scale_factor_1, scale_factor_2, scale_factor_3, pitch_kick, mleft_1, mleft_2, mleft_3, count, delta_v_total, delta_v_circularization

%-------- Variables for Inputs for Function defined here (in kg,m,s format)

One_dv_params = readtable('LVMasses\OnedVParameters.csv'); % load Mission 1 dv requirements csv into table
Two_dv_params = readtable('LVMasses\TwodVParameters.csv'); % load Mission 2 dv requirements csv into table 
Latona1_params = readtable('LVTrajectory\Latona-1TrajReqs.csv');
Latona2_params = readtable('LVTrajectory\Latona-2TrajReqs.csv');
Minerva1_params = readtable('LVTrajectory\Minerva-1TrajReqs.csv');
Minerva2_params = readtable('LVTrajectory\Minerva-2TrajReqs.csv');
Zephyr1_params = readtable('LVTrajectory\Zephyr-1TrajReqs.csv');
Zephyr2_params = readtable('LVTrajectory\Zephyr-2TrajReqs.csv');
%--Zephyr Inputs--
grav_sea_level = 9.80665;
Mass_Payload_1 = 30;
Mass_Payload_2 = 95;

% Zephyr_1_Cd2 = 0.2;
% Zephyr_1_Radius2 = 0.45; 
% Zephyr_1_mi2 = 1336.692; 
% Zephyr_1_mf2 = 143.7361 + Mass_Payload_1; 
% Zephyr_1_mprop2 = Zephyr_1_mi2 - Zephyr_1_mf2;
% Zephyr_1_T2 = 1.05 * (Zephyr_1_mi2) * grav_sea_level;
% Zephyr_1_Isp2 = 359.1;
% 
% Zephyr_1_Cd1 = 0.2;
% Zephyr_1_Radius1 = 0.45;
% Zephyr_1_mi1 = 4338.224; 
% Zephyr_1_mf1 = 447.2046; 
% Zephyr_1_mprop1 = Zephyr_1_mi1 - Zephyr_1_mf1;
% Zephyr_1_T1 = 1.4 * (Zephyr_1_mi1 + Zephyr_1_mi2) * grav_sea_level;
% Zephyr_1_Isp1 = 296.1;
% 
% %---
% 
% Zephyr_2_Cd3 = 0.2;
% Zephyr_2_Radius3 = 0.42;
% Zephyr_2_mi3 = 426.4761;
% Zephyr_2_mf3 = 36.46237 + Mass_Payload_2; 
% Zephyr_2_mprop3 = Zephyr_2_mi3 - Zephyr_2_mf3;
% Zephyr_2_T3 = 0.9 * (Zephyr_2_mi3) * grav_sea_level;
% Zephyr_2_Isp3 = 359.1;
% 
% Zephyr_2_Cd2 = 0.2;
% Zephyr_2_Radius2 = 0.45;
% Zephyr_2_mi2 = 1306.692;
% Zephyr_2_mf2 = 143.7361;
% Zephyr_2_mprop2 = Zephyr_2_mi2 - Zephyr_2_mf2;
% Zephyr_2_T2 = 1.05 * (Zephyr_2_mi2 + Zephyr_2_mi3) * grav_sea_level;
% Zephyr_2_Isp2 = 359.1;
% 
% Zephyr_2_Cd1 = 0.2;
% Zephyr_2_Radius1 = 0.45;
% Zephyr_2_mi1 = 4338.224; 
% Zephyr_2_mf1 = 477.2046; 
% Zephyr_2_mprop1 = Zephyr_2_mi1 - Zephyr_2_mf1;
% Zephyr_2_T1 = 1.4 * (Zephyr_2_mi1 + Zephyr_2_mi2 + Zephyr_2_mi3) * grav_sea_level;
% Zephyr_2_Isp1 = 296.1;

Zephyr_1_Cd2 = 0.2;
Zephyr_1_Radius2 = Zephyr1_params{2,2}; 
Zephyr_1_mi2 = Zephyr1_params{2,3};
Zephyr_1_mf2 = Zephyr1_params{2,4} + Mass_Payload_1;
Zephyr_1_mprop2 = Zephyr_1_mi2 - Zephyr_1_mf2;
Zephyr_1_T2 = 1.05 * (Zephyr_1_mi2) * grav_sea_level;
Zephyr_1_Isp2 = Zephyr1_params{2,8};

Zephyr_1_Cd1 = 0.2;
Zephyr_1_Radius1 = Zephyr1_params{1,2};
Zephyr_1_mi1 = Zephyr1_params{1,3};
Zephyr_1_mf1 = Zephyr1_params{1,4};
Zephyr_1_mprop1 = Zephyr_1_mi1 - Zephyr_1_mf1;
Zephyr_1_T1 = 1.4 * (Zephyr_1_mi1 + Zephyr_1_mi2) * grav_sea_level;
Zephyr_1_Isp1 = Zephyr1_params{1,8};

%---

Zephyr_2_Cd3 = 0.2;
Zephyr_2_Radius3 = Zephyr2_params{3,2};
Zephyr_2_mi3 = Zephyr2_params{3,3};
Zephyr_2_mf3 = Zephyr2_params{3,4} + Mass_Payload_2;
Zephyr_2_mprop3 = Zephyr_2_mi3 - Zephyr_2_mf3;
Zephyr_2_T3 = 0.9 * (Zephyr_2_mi3) * grav_sea_level;
Zephyr_2_Isp3 = Zephyr2_params{3,8};

Zephyr_2_Cd2 = 0.2;
Zephyr_2_Radius2 = Zephyr2_params{2,2};
Zephyr_2_mi2 = Zephyr2_params{2,3};
%disp(Zephyr2_params{2,3});
Zephyr_2_mf2 = Zephyr2_params{2,4};
Zephyr_2_mprop2 = Zephyr_2_mi2 - Zephyr_2_mf2;
Zephyr_2_T2 = 1.05 * (Zephyr_2_mi2 + Zephyr_2_mi3) * grav_sea_level;
Zephyr_2_Isp2 = Zephyr2_params{2,8};

Zephyr_2_Cd1 = 0.2;
Zephyr_2_Radius1 = Zephyr2_params{1,2};
Zephyr_2_mi1 = Zephyr2_params{1,3};
%disp(Zephyr2_params{1,3});
Zephyr_2_mf1 = Zephyr2_params{1,4};
Zephyr_2_mprop1 = Zephyr_2_mi1 - Zephyr_2_mf1;
Zephyr_2_T1 = 1.4 * (Zephyr_2_mi1 + Zephyr_2_mi2 + Zephyr_2_mi3) * grav_sea_level;
Zephyr_2_Isp1 = Zephyr2_params{1,8};

%----Using Trajectory Function with inputvariables above

Zephyr_1_Trajectories = Trajectory_TL_0407('Zephyr',[Zephyr_1_Cd1; Zephyr_1_Radius1; Zephyr_1_mi1; Zephyr_1_mf1; Zephyr_1_T1 ; Zephyr_1_Isp1] , [Zephyr_1_Cd2; Zephyr_1_Radius2; Zephyr_1_mi2; Zephyr_1_mf2; Zephyr_1_T2; Zephyr_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, 28.50);
Zephyr_2_Trajectories = Trajectory_TL_0407('Zephyr',[Zephyr_2_Cd1; Zephyr_2_Radius1; Zephyr_2_mi1; Zephyr_2_mf1; Zephyr_2_T1 ; Zephyr_2_Isp1] , [Zephyr_2_Cd2; Zephyr_2_Radius2; Zephyr_2_mi2; Zephyr_2_mf2; Zephyr_2_T2; Zephyr_2_Isp2] , [Zephyr_2_Cd3; Zephyr_2_Radius3; Zephyr_2_mi3; Zephyr_2_mf3; Zephyr_2_T3; Zephyr_2_Isp3] , 2, 34.60);

[Minimum_Values_1 , Index_1] = min(Zephyr_1_Trajectories, [], 1);
[Minimum_Values_2 , Index_2] = min(Zephyr_2_Trajectories, [], 1);

Optimal_Zephyr_1_Result = Zephyr_1_Trajectories(Index_1(10),:);
Optimal_Zephyr_2_Result = Zephyr_2_Trajectories(Index_2(10),:);

%-------------------------------------------------------------------------
%This code runs the Trajectory function once for when the optimal
%trajectory is found and graphs it

Optimal_Zephyr_1_Trajectory = Optimal_Trajectory('Zephyr',[Zephyr_1_Cd1; Zephyr_1_Radius1; Zephyr_1_mi1; Zephyr_1_mf1; Zephyr_1_T1 ; Zephyr_1_Isp1] , [Zephyr_1_Cd2; Zephyr_1_Radius2; Zephyr_1_mi2; Zephyr_1_mf2; Zephyr_1_T2; Zephyr_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, 28.50,Optimal_Zephyr_1_Result);
Optimal_Zephyr_2_Trajectory = Optimal_Trajectory('Zephyr',[Zephyr_2_Cd1; Zephyr_2_Radius1; Zephyr_2_mi1; Zephyr_2_mf1; Zephyr_2_T1 ; Zephyr_2_Isp1] , [Zephyr_2_Cd2; Zephyr_2_Radius2; Zephyr_2_mi2; Zephyr_2_mf2; Zephyr_2_T2; Zephyr_2_Isp2] , [Zephyr_2_Cd3; Zephyr_2_Radius3; Zephyr_2_mi3; Zephyr_2_mf3; Zephyr_2_T3; Zephyr_2_Isp3] , 2, 34.60,Optimal_Zephyr_2_Result);

Optimal_Zephyr_1_Trajectory = Optimal_Zephyr_1_Trajectory(1:Minimum_Values_1(8), :);
Optimal_Zephyr_2_Trajectory = Optimal_Zephyr_2_Trajectory(1:Minimum_Values_2(8), :);

scales1 = Optimal_Zephyr_1_Result(1:7);
scales2 = Optimal_Zephyr_2_Result(1:7);

plotTrajZ('Zephyr',[Zephyr_1_Cd1; Zephyr_1_Radius1; Zephyr_1_mi1; Zephyr_1_mf1; Zephyr_1_T1 ; Zephyr_1_Isp1] , [Zephyr_1_Cd2; Zephyr_1_Radius2; Zephyr_1_mi2; Zephyr_1_mf2; Zephyr_1_T2; Zephyr_1_Isp2] , [0; 0; 0; 0; 0; 0], 1, scales1);
plotTrajZ('Zephyr',[Zephyr_2_Cd1; Zephyr_2_Radius1; Zephyr_2_mi1; Zephyr_2_mf1; Zephyr_2_T1 ; Zephyr_2_Isp1] , [Zephyr_2_Cd2; Zephyr_2_Radius2; Zephyr_2_mi2; Zephyr_2_mf2; Zephyr_2_T2; Zephyr_2_Isp2] , [Zephyr_2_Cd3; Zephyr_2_Radius3; Zephyr_2_mi3; Zephyr_2_mf3; Zephyr_2_T3; Zephyr_2_Isp3] , 2, scales2);

save('Trajectory_TL_0407_Workspace.mat');
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
                                stage1_S = pi*stage1_Radius^2; %m^2                      %Stage 1 Cross-Sectional Area
                                stage1_Thrust = scale_factor_1*step1(5); %N              %Stage 1 Thrust
                                stage1_Isp = step1(6); %s                                %Stage 1 Isp
                                stage1_mdot = stage1_Thrust/(stage1_Isp*g0); %kg/s       %Stage 1 Mass Flow
                                stage1_mf = step1(4)+step2(3)+step3(3); %kg              %Stage 1 Final Mass (Stage 1 Struct, Stage 2 Fuel, Stage 2 Struct)     
                                
                                
                                stage2_Cd = step2(1);                                    %Stage 2 Coefficient of Drag
                                stage2_Radius = step2(2); %m                             %Stage 2 Radius
                                stage2_mi = step2(3)+step3(3); %kg                       %Stage 2 Initial Mass (Stage l Struct, Stage 1 Fuel, Stage 2 Struct, Stage 2 Fuel)
                                stage2_mf = step2(4)+step3(3); %kg                       %Stage 2 Final Mass (Stage 1 Struct, Stage 2 Fuel, Stage 2 Struct)
                                stage2_S = pi*stage2_Radius^2; %m^2                      %Stage 2 Cross-Sectional Area
                                stage2_Thrust = scale_factor_2*step2(5); %N                          %Stage 2 Thrust
                                stage2_Isp = step2(6); %s                                %Stage 2 Isp
                                stage2_mdot = stage2_Thrust/(stage2_Isp*g0); %kg/s       %Stage 2 Mass Flow
                                
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
                                        m(count+1) = m(count) - stage1_mdot*dt; 
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

                                %% Hohmann Transfer
                                r_periapsis = h(count) + R_earth;
                                r_apoapsis = rf;
                                eccentricity_transfer = (r_apoapsis - r_periapsis)/(r_apoapsis + r_periapsis);

                                if inc < 90
                                    delta_v_circularization = sqrt(mu/r_periapsis) - (v(count) + v_ls);
                                else
                                    delta_v_circularization = sqrt(mu/r_periapsis) - (v(count) - v_ls);
                                end
                                delta_v_1 = sqrt(mu/r_periapsis) * (sqrt((2*r_apoapsis)/(r_periapsis+r_apoapsis))-1);
                                delta_v_2 = sqrt(mu/r_apoapsis) * (1-sqrt((2*r_periapsis)/(r_periapsis+r_apoapsis)));

                                delta_v_total = delta_v_1 + delta_v_2 + delta_v_circularization;

                                Optimal_Result = [t,v,m,h,x,gamma,a,q];
    end


