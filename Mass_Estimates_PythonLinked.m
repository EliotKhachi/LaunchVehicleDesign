clc
clear
%==========================================================================
% This code develops the initial Mass Estimations
% 
% Developed by: Ape Lincoln
%
% Editors notes: 
% Change only the values in the Inputs category
% First values in each array is for Minerva
% Second Value in each array is for Latona
%
%==========================================================================
% known variables
g0 = 9.80665;                             % m/s^2
m_payload = [30 95];                      % Payload for mission 1 and 2 respectively (kg)

%==========================================================================
% input
% 1 row, 7 columns
One_dv_params = readtable('LVMasses\OnedVParameters.csv'); % load Mission 1 dv requirements csv into table
Two_dv_params = readtable('LVMasses\TwodVParameters.csv'); % load Mission 2 dv requirements csv into table 
% col 1 - dv needed/ideal, col 2 - dv design, col 3 - dv plane,
% col 4 - dv grav loss, col 5 -dv drag loss, col 6 - apo kick,
% col 7 - dv manuevers

%V_loss = [1759.55, 1759.55];               % Velocity from losses for Minerva and Latona (m/s)
V_loss = [One_dv_params{1,4}+One_dv_params{1,5}, Two_dv_params{1,4}+Two_dv_params{1,5}];         % Velocity from losses for Minerva and Latona (m/s)
%V_plane = 1304.936;                                % Mission 2 plane change DV
V_plane = Two_dv_params{1,3}; 
%V_Orbital = [85.34026, 98.92758+V_plane];          % Other DV losses for LV
V_Orbital = [One_dv_params{1,6}, Two_dv_params{1,6}+V_plane];  % Other DV losses for LV
%V_Recover = 0.0000;  
V_Recover = Two_dv_params{1,7};                     % Recovery DV
%V_ideal = [7687.89945 7829.83702];                 % Ideal velocity of mission 1 and 2 Respectively(m/s)
V_ideal = [One_dv_params{1,1}, Two_dv_params{1,1}]; % Ideal velocity of mission 1 and 2 Respectively(m/s)
% sigma_1 = [0.11 0.08];                    % Stage 1 structural mass fraction
% sigma_2 = [0.11 0.11];                    % Stage 2 structural mass fraction
% sigma_3 = [0.11 0.08];                    % Stage 3/0 structural mass fraction
sigma_1 = [0.11 0.13];                    % Stage 1 structural mass fraction
sigma_2 = [0.11 0.2];                    % Stage 2 structural mass fraction
sigma_3 = [0.11 0.13];                    % Stage 3/0 structural mass fraction
isp_1 = [296.1 265];                      % Stage 1 Isp (s)
isp_2 = [359.1 380];                      % Stage 2 Isp (s)
isp_3 = [359.1 265];                      % Stage 3 Isp (s)

%==========================================================================
% code
for i=1:2                                 % LV Selection
    for j=1:2                             % Mission Selection
        if i==1
            if j==1
                m_0_1 = [];
                DV_split = [];
                Minerva_1DV = [];
                Minerva_2DV = [];
                for DV_percent1 = 0.15:0.05:0.50
                for DV_percent3 = 0.15:0.05:0.50
                    DV_percent2 = 1-(DV_percent1+DV_percent3);
                    
                    V_stage3 = V_ideal(j+1)*DV_percent3;
                    DV_req3 = V_stage3+V_Orbital(j+1);
                    mu_3 = exp((DV_req3)/(isp_3(i)*g0));
                    m_p3 = m_payload(j+1)*(((mu_3-1)*(1-sigma_3(i)))./(1-(sigma_3(i)*mu_3)));
                    m_s3 = m_p3*((sigma_3(i))/(1-sigma_3(i)));
                    m_03 = m_p3+m_s3+m_payload(j+1);
            
                    V_stage2 = V_ideal(j+1)*DV_percent2;
                    DV_req2 = V_stage2;
                    mu_2 = exp((DV_req2)/(isp_2(i)*g0));
                    m_p2 = m_03.*(((mu_2-1)*(1-sigma_2(i)))./(1-(sigma_2(i)*mu_2)));
                    m_s2 = m_p2*((sigma_2(i))/(1-sigma_2(i)));
                    m_02 = m_p2+m_s2;
            
                    V_stage1 = V_ideal(j+1)*DV_percent1;
                    DV_req1 = V_stage1+V_Recover+V_loss(j+1);
                    mu_1 = exp((DV_req1)/(isp_1(i)*g0));
                    m_p1 = (m_02+m_03).*(((mu_1-1)*(1-sigma_1(i)))./(1-(sigma_1(i)*mu_1)));
                    m_s1 = m_p1*((sigma_1(i))/(1-sigma_1(i)));
                    m_01 = m_p1+m_s1;
                    
                    m_1i = m_01+m_02+m_payload(j);
                    m_1f = m_s1+m_02+m_payload(j);
                    DV_req1_1 = isp_1(i)*g0*log(m_1i/m_1f);
                    m_2i = m_02+m_payload(j);
                    m_2f = m_s2+m_payload(j);
                    DV_req2_1 = isp_2(i)*g0*log(m_2i/m_2f);
                    DV_req_1 = DV_req2_1+DV_req1_1;
                    if DV_req_1 >= V_ideal(j)+V_loss(j)+V_Recover+V_Orbital(j)
                        m_0_1 = [m_0_1,m_1i];
                        DV_split = [DV_split,[DV_percent3;DV_percent2;DV_percent1]];
                    end
                    
                end
                end
                [m_optimal,I] = min(m_0_1);
                Minerva_mass = m_optimal;
                Minerva_split = DV_split(:,I);
                
                DV_percent3 = Minerva_split(1);
                DV_percent2 = Minerva_split(2);
                DV_percent1 = Minerva_split(3);
                
                V_stage3 = V_ideal(j+1)*DV_percent3;
                DV_req3 = V_stage3+V_Orbital(j+1);
                mu_3 = exp((DV_req3)/(isp_3(i)*g0));
                m_p3 = m_payload(j+1)*(((mu_3-1)*(1-sigma_3(i)))./(1-(sigma_3(i)*mu_3)));
                m_s3 = m_p3*((sigma_3(i))/(1-sigma_3(i)));
                m_03 = m_p3+m_s3+m_payload(j+1);
        
                V_stage2 = V_ideal(j+1)*DV_percent2;
                DV_req2 = V_stage2;
                mu_2 = exp((DV_req2)/(isp_2(i)*g0));
                m_p2 = m_03.*(((mu_2-1)*(1-sigma_2(i)))./(1-(sigma_2(i)*mu_2)));
                m_s2 = m_p2*((sigma_2(i))/(1-sigma_2(i)));
                m_02 = m_p2+m_s2;
        
                V_stage1 = V_ideal(j+1)*DV_percent1;
                DV_req1 = V_stage1+V_Recover+V_loss(j+1);
                mu_1 = exp((DV_req1)/(isp_1(i)*g0));
                m_p1 = (m_02+m_03).*(((mu_1-1)*(1-sigma_1(i)))./(1-(sigma_1(i)*mu_1)));
                m_s1 = m_p1*((sigma_1(i))/(1-sigma_1(i)));
                m_01 = m_p1+m_s1;
                
                m_0 = m_01+m_02+m_03;
                
                Minerva_2DV = {'Minerva Mission 2','Stage 3','Stage 2','Stage 1';'Sigma',sigma_3(i),sigma_2(i),sigma_1(i);'Isp',isp_3(i),isp_2(i),isp_1(i);'%DV Split',DV_percent3,DV_percent2,DV_percent1;'Propellant Mass',m_p3,m_p2,m_p1;'Structural mass',m_s3,m_s2,m_s1;'Total Stage Mass',m_03,m_02,m_01;'Total Mass',' ',' ',m_0};
                Minerva2 = [sigma_3(i),sigma_2(i),sigma_1(i);isp_3(i),isp_2(i),isp_1(i);DV_percent3,DV_percent2,DV_percent1;m_p3,m_p2,m_p1;m_s3,m_s2,m_s1;m_03,m_02,m_01;0,0,m_0];
                writecell(Minerva_2DV,'LVMasses\Minerva-2MassEstimate.csv')
                
                m_p3 = 0;
                m_s3 = 0;
                m_03 = 0;
                DV_percent3 = 0;
                
                mu_2 = (m_p2+(m_payload(j)*(1-sigma_2(i))))/((m_p2*sigma_2(i))+(m_payload(j)*(1-sigma_2(i))));
                DV_req2 = log(mu_2)*isp_2(i)*g0;
                V_stage2 = DV_req2-V_Orbital(j);
                DV_percent2 = V_stage2/V_ideal(j);
                m_02 = m_02+m_payload(j);
                
                mu_1 = (m_p1+(m_02*(1-sigma_1(i))))/((m_p1*sigma_1(i))+(m_02*(1-sigma_1(i))));
                DV_req1 = log(mu_1)*isp_1(i)*g0;
                V_stage1 = DV_req1-V_Recover-V_loss(j);
                DV_percent1 = V_stage1/V_ideal(j);
                
                m_0 = m_01+m_02+m_payload(j);
                
                Minerva_1DV = {'Minerva Mission 1','Stage 3','Stage 2','Stage 1';'Sigma',sigma_3(i),sigma_2(i),sigma_1(i);'Isp',isp_3(i),isp_2(i),isp_1(i);'%DV Split',DV_percent3,DV_percent2,DV_percent1;'Propellant Mass',m_p3,m_p2,m_p1;'Structural mass',m_s3,m_s2,m_s1;'Total Stage Mass',m_03,m_02,m_01;'Total Mass',' ',' ',m_0};
                Minerva1 = [sigma_3(i),sigma_2(i),sigma_1(i);isp_3(i),isp_2(i),isp_1(i);DV_percent3,DV_percent2,DV_percent1;m_p3,m_p2,m_p1;m_s3,m_s2,m_s1;m_03,m_02,m_01;0,0,m_0];
                writecell(Minerva_1DV,'LVMasses\Minerva-1MassEstimate.csv')
            end
        else
        if i==2
            if j==1
                DV_percent1 = 0.15:0.05:0.95;
                DV_percent2 = 1-DV_percent1;
                DV_percent0 = 1-(DV_percent1+DV_percent2);
            
                V_stage2 = V_ideal(j)*DV_percent2;
                DV_req2 = V_stage2+V_Orbital(j);
                mu_2 = exp((DV_req2)/(isp_2(i)*g0));
                m_p2 = m_payload(j)*(((mu_2-1)*(1-sigma_2(i)))./(1-(sigma_2(i)*mu_2)));
                m_s2 = m_p2*((sigma_2(i))/(1-sigma_2(i)));
                m_02 = m_p2+m_s2+m_payload(j);
            
                V_stage1 = V_ideal(j)*DV_percent1;
                DV_req1 = V_stage1+V_loss(i);
                mu_1 = exp((DV_req1)/(isp_1(i)*g0));
                m_p1 = m_02.*(((mu_1-1)*(1-sigma_1(i)))./(1-(sigma_1(i)*mu_1)));
                m_s1 = m_p1*((sigma_1(i))/(1-sigma_1(i)));
                m_01 = m_p1+m_s1;
            
                V_stage0 = 0;
                DV_req0 = 0;
                mu_0 = 0;
                m_p0 = 0;
                m_s0 = 0;
                m_00 = 0;
                
                m_0 = m_02+m_01+m_00;
                for k=1:length(DV_percent1)
                   if m_0(k)<0
                       m_0(k)=NaN;
                   end
                end
                [m_optimal,I] = min(m_0,[],'omitnan');
            
                Latona_1DV = {'Latona Mission 1','Stage 2','Stage 1','Stage 0';'Sigma',sigma_2(i),sigma_1(i),sigma_3(i);'Isp',isp_2(i),isp_1(i),isp_3(i);'%DV Split',DV_percent2(I),DV_percent1(I),DV_percent0(I);'Propellant Mass',m_p2(I),m_p1(I),m_p0;'Structural mass',m_s2(I),m_s1(I),m_s0;'Total Stage Mass',m_02(I),m_01(I),m_00;'Total Mass',' ',' ',m_0(I)};
                Latona1 = [sigma_2(i),sigma_1(i),sigma_3(i);isp_2(i),isp_1(i),isp_3(i);DV_percent2(I),DV_percent1(I),DV_percent0(I);m_p2(I),m_p1(I),m_p0;m_s2(I),m_s1(I),m_s0;m_02(I),m_01(I),m_00;0,0,m_0(I)];
                writecell(Latona_1DV,'LVMasses\Latona-1MassEstimate.csv')
            else
            if j==2
                syms DV_p1 DV_p0
                m_p2 = m_p2(I);
                m_s2 = m_s2(I);
                m_02 = m_p2+m_s2+m_payload(j);
                mu_2 = (m_p2+(m_payload(j)*(1-sigma_2(i))))/((sigma_2(i)*m_p2)+(m_payload(j)*(1-sigma_2(i))));
                DV_req2 = log(mu_2)*isp_2(i)*g0;
                V_stage2 = DV_req2-V_Orbital(j);
                DV_percent2 = V_stage2/V_ideal(j);
            
                m_p1 = m_p1(I);
                m_s1 = m_s1(I);
                m_01 = m_p1+m_s1;
                mu_1 = (m_p1+(m_02*(1-sigma_1(i))))/((sigma_1(i)*m_p1)+(m_02*(1-sigma_1(i))));
                DV_req1 = log(mu_1)*isp_1(i)*g0;
                V_losspercent1 = (DV_p1)/(DV_p1+DV_p0);
                V_stage1 = DV_req1-(V_losspercent1*V_loss(j));
                eqn1 = DV_p1 == V_stage1/V_ideal(j);
            
                eqn2 = DV_p0 == 1-(DV_p1+DV_percent2);
                s = vpasolve([eqn1 eqn2],[DV_p1 DV_p0],[0 inf]);
                DV_percent1 = double(s.DV_p1);
                DV_percent0 = double(s.DV_p0);
                
                V_losspercent1 = (DV_percent1)/(DV_percent1+DV_percent0);
                V_stage1 = DV_req1-(V_losspercent1*V_loss(j));
                
                r_b = 1;
                V_stage0 = V_ideal(j)*DV_percent0;
                V_losspercent0 = (DV_percent0)/(DV_percent1+DV_percent0);
                DV_req0 = V_stage0+(V_losspercent0*DV_percent0);
                mu_0 = exp((DV_req0)/(isp_3(i)*g0));
                m_p0 = (m_01+m_02)*(((mu_0-1)*(1-sigma_3(i)))./(1-(sigma_3(i)*mu_0)));
                m_s0 = m_p0*((sigma_3(i))/(1-sigma_3(i)));
                m_00 = m_p0+m_s0;
                
                m_0 = m_02+m_01+m_00;
                
                TW = 1.4;
                T = TW*m_0*g0;
                mdot = T/g0;
                t_b = m_p0/mdot;
                
                V_idealtotal = V_stage2+V_stage1+V_stage0;
                V_losstotal = V_loss(j)*(V_losspercent0+V_losspercent1);
                
                Latona_2DV = {'Latona Mission 2','Stage 2','Stage 1','Stage 0';'Sigma',sigma_2(i),sigma_1(i),sigma_3(i);'Isp',isp_2(i),isp_1(i),isp_3(i);'%DV Split',DV_percent2,DV_percent1,DV_percent0;'Propellant Mass',m_p2,m_p1,m_p0;'Structural mass',m_s2,m_s1,m_s0;'Total Stage Mass',m_02,m_01,m_00;'Total Mass',' ',' ',m_0};
                Latona2 = [sigma_2(i),sigma_1(i),sigma_3(i);isp_2(i),isp_1(i),isp_3(i);DV_percent2,DV_percent1,DV_percent0;m_p2,m_p1,m_p0;m_s2,m_s1,m_s0;m_02,m_01,m_00;0,0,m_0];
                writecell(Latona_2DV,'LVMasses\Latona-2MassEstimate.csv')
            end
            end
        end
        end
    end
end