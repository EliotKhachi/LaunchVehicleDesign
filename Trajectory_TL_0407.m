%------------Editor's Notes, Updates, and Changes ------------------------
%Changed to focus on Zephyr

%----------- How to Use --------------------------------------------------
%name is a string with the name of the rocket being explore: Zephyr

%step1 is a column array containg the Cd, Radius, mass initial (fuel + struct), mass
%final (struct), thrust, and Isp of the first stage (in metric)

%step2 is a column array containg the Cd, Radius, mass initial (fuel + struct), mass
%final (struct), thrust, and Isp of the second stage (in metric)

%step is a column array containg the Cd, Radius, mass initial (fuel + struct), mass
%final (struct), thrust, and Isp of the third stage (in metric)

%mission is a double which desribes the mission we're testing: 1 or 2

%launch_latitude is a double which describes the latitude of the launch 
%site (in degrees)

%Example function inputs: Trajectory_<name>_<date>('Zephyr',[0.2; 0.8; 1500; 1000; 30000; 280],[0.2; 0.8; 750; 100; 25000; 275], [0.0; 0.0; 0; 0; 0; 0],1,35)
%Example function inputs based off HW7: Trajectory_<name>_<date>('Zephyr',[0.2; 5; 3000000; 100000; 33000000; 450],[0; 0; 0; 0; 0; 0], [0.0; 0.0; 0; 0; 0; 0], 1, 28.45)
%Note: The above inputs will yeild erronous output. Only use to test
%functionality of the function.

%Vandenberg: Lat = 34.60 deg.
%Kennedy Space Center: Lat = 28.50 deg.

function Results = Trajectory_TL_0129(name, step1, step2, step3, mission, launch_latitude)
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
    stage1_Thrust_Scaling_Factor = .8:.05:1.2;
    %stage1_Thrust_Scaling_Factor = 0.9;
    
    stage2_Thrust_Scaling_Factor = .7:.05:1.1;
    %stage2_Thrust_Scaling_Factor = 0.6;
    
    stage3_Thrust_Scaling_Factor = .7:.05:1.1;
    %stage3_Thrust_Scaling_Factor = 0.5;
    
    pitchkick_spread = 0.05:.05:1.5;
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
    q0 = 0; %Pa                                         %Initial Dynamic Pressure
    Results = [];
    Record = [];
    %% Simulation
    
    for scale_factor_1 = stage1_Thrust_Scaling_Factor
        for scale_factor_2 = stage2_Thrust_Scaling_Factor
            for scale_factor_3 = stage3_Thrust_Scaling_Factor
                for pitch_kick = pitchkick_spread
                    for mleft_1 = stage1_mleft
                        for mleft_2 = stage2_mleft
                            for mleft_3 = stage3_mleft
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
                                stage3_mdot = 0; %kg/s                                   %Stage 3 Mass Flow
                                
                                
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

                                if stage3_mi == 0
                                    MR_H_Transfer = exp(delta_v_total / (g0*stage2_Isp));
                                    if MR_H_Transfer > ((stage2_mf+mleft_2*(stage2_mi-stage2_mf))/(stage2_mf))
                                        check = 0;
                                    else
                                        check = 1;
                                    end
                                else
                                    MR_H_Transfer = exp(delta_v_total / (g0*stage3_Isp));
                                    if MR_H_Transfer > ((stage3_mf+mleft_3*(stage3_mi-stage3_mf))/(stage3_mf))
                                        check = 0;
                                    else
                                        check = 1;
                                    end
                                end
                                
                                if mission == 2
                                   if h(count) >= 450000 || h(count) < 0
                                   check = 0; 
                                   end    
                                else
                                if h(count) >= 250000 || h(count) < 0
                                    check = 0; 
                                end
                                end

                                if gamma > deg2rad(1)
                                      check = 0;
                                end

                                if delta_v_total <= 0 || delta_v_circularization <= 0 
                                    check = 0;
                                end
                                
                                if delta_v_circularization >= 500 && strcmp(name, 'Zephyr') && mission == 1
                                   check = 0; 
                                end
                                
                                if delta_v_circularization >= 1000 && strcmp(name, 'Zephyr') && mission == 2
                                   check = 0; 
                                end 
                                
                                v_circ_final_2 = sqrt( mu / ( R_earth + 550000 ));
                                inc_change_2 = 2 * v_circ_final_2 * sin(deg2rad(10)/2);
                                dv_check = g0 * stage3_Isp * log((mleft_3*(stage3_mi-stage3_mf) + stage3_mf) / (stage3_mf));
                                
                                if dv_check <= inc_change_2 && mission == 2
                                   check = 0;
                                end
                                
                                if check == 1
                                    Results = [Results; scale_factor_1, scale_factor_2, scale_factor_3, pitch_kick, mleft_1, mleft_2, mleft_3, count, delta_v_total, delta_v_circularization,h(count),v(count)];
                                end
                                
                                Record = [Record; scale_factor_1, scale_factor_2, scale_factor_3, pitch_kick, mleft_1, mleft_2, mleft_3, count, delta_v_total, delta_v_circularization,h(count),v(count)];

                            end
                        end
                    end
                end
            end
        end
    end
end


