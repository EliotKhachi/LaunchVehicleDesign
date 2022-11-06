function plotTraj(name, step1, step2, step3, mission, scales)
%% Earth and Launch Site Inputs
    g0 = 9.80665; %m/s^2                Local Gravity at S.L.
    R_earth = 6378000; %m               Radius of Earth
    h0 = 7640; %m                       Scale Height (Wikipedia: https://en.wikipedia.org/wiki/Scale_height)
    rho0 = 1.225; %kg/m^3               Air Density at S.L.
    
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
    %% Simulation
    for scale_factor_1 = scales(1)
        for scale_factor_2 = scales(2)
            for scale_factor_3 = scales(3)
                for pitch_kick = scales(4)
                    for mleft_1 = scales(5)
                        for mleft_2 = scales(6)
                            for mleft_3 = scales(7)
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
                                    stage1_Booster_mass_indv(1) = step1(3)/4;
                                    stage1_Booster_mdot_indv = stage1_mdot/4;
                                    stage1_Main_mass = step2(3);
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
                                %Circular Orbit Velocity of 150,000 the
                                %minimum orbit altitude is 7815 m/s
                                count = 1;                                               %Begin Count Here
                                while gamma(count) > deg2rad(1)  
                                    t(count+1) = t(count)+1;
                                    if m(count) > stage1_mf+(mleft_1*(stage1_mi-stage1_mf))%The Order of calculating these terms are important
                                        v(count+1) = v(count) + a(count)*dt;
                                        if mission == 2 && strcmp(name, 'Latona')
                                           m(count+1) = m(count) - stage1_mdot*dt - stage2_mdot*dt;
                                           stage1_Main_mass(count+1) = stage1_Main_mass(count) - stage2_mdot*dt;
                                           stage1_Booster_mass_indv(count+1) = stage1_Booster_mass_indv(count) - stage1_Booster_mdot_indv*dt;
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
                                yyaxis left
                                plot(t(1:count), h(1:count)/1000, t(1:count), rad2deg(gamma(1:count)), t(1:count), a(1:count), '-.')
                                xlabel('Time (s)')
                                ylabel('Altitude (km), Flight Path Angle (Deg), Acceleration (m/s^2)')
                                ylim([0 inf])

                                yyaxis right
                                plot(t(1:count), v(1:count), t(1:count), q(1:count)/10)
                                ylabel('Velocity (m/s), dynamic pressure (Pa/10)')
                                legend('altitude', 'gamma', 'acceleration', 'velocity', 'dynamic pressure', 'Location', 'northwest')
                                title(strcat(name, {' '}, string(mission), ' Trajectory'))
                                ylim([0 inf])
                                grid on
                                
                                [~, ind] = max(q);
                                xline(t(ind))
                                
                                if mission == 2 && strcmp(name, 'Latona')
                                    T = {'Time (s)', 'Thrust (N)', 'Max-q (Pa)', 'Velocity (m/s)', 'Mass Burned (kg)', 'Stap-on Booster Masses (kg)', 'Main Booster Mass (kg)', 'Height (m)', 'Gamma (rad)', 'Air Density (kg/m^3)';...
                                     t(ind), stage1_Thrust, q(ind), v(ind), m(1) - m(ind), stage1_Booster_mass_indv(ind), stage1_Main_mass(ind), h(ind), gamma(ind), rho(ind)};
                                 
                                else
                                    T = {'Time (s)', 'Thrust (N)', 'Max-q (Pa)', 'Velocity (m/s)', 'Mass Burned (kg)', 'Height (m)', 'Gamma (rad)', 'Air Density (kg/m^3)';...
                                     t(ind), stage1_Thrust, q(ind), v(ind), m(1) - m(ind), h(ind), gamma(ind), rho(ind)};
                                 
                                end
                                 
                                 
                                writecell(T, "LVMasses\Max Q Conditions_" + string(name) + "-" + mission + ".csv")
                                name = name + " " + string(mission);
                                saveas(figure(1), strcat(name, '.png'))
                            end
                        end
                    end
                end
            end
        end
    end
end