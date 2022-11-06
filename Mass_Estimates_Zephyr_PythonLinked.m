clc
clear
%==========================================================================
% This code develops the initial Mass Estimations
% 
% Developed by: Ape Lincoln
%
% Editors notes: 
% Change only the values in the Inputs category
% First values in each array is for Zephyr
% Second Value in each array is for Latona
%
%==========================================================================
% known variables
g0 = 9.80665;                             % m/s^2
m_payload = [30 95];                      % Payload for mission 1 and 2 respectively (kg)

%==========================================================================
% input
One_dv_params = readtable('LVMasses\OnedVParameters.csv'); % load Mission 1 dv requirements csv into table
Two_dv_params = readtable('LVMasses\TwodVParameters.csv'); % load Mission 2 dv requirements csv into table 

% V_loss = [1759.55048 1759.55048];         % Velocity from losses for Zephyr and Latona (m/s)
% V_plane = 1304.93618;                     % Mission 2 plane change DV
% V_Orbital = [85.34026 98.92758+V_plane];  % Other DV losses for LV
% V_Recover = 0.00000;                      % Recovery DV
% V_ideal = [7687.89945 7829.83702];        % Ideal velocity of mission 1 and 2 Respectively(m/s)

V_loss = [One_dv_params{1,4}+One_dv_params{1,5}, Two_dv_params{1,4}+Two_dv_params{1,5}];         % Velocity from losses for Minerva and Latona (m/s)
V_plane = Two_dv_params{1,3}; 
V_Orbital = [One_dv_params{1,6}, Two_dv_params{1,6}+V_plane];  % Other DV losses for LV
V_Recover = Two_dv_params{1,7};                     % Recovery DV
V_ideal = [One_dv_params{1,1}, Two_dv_params{1,1}]; % Ideal velocity of mission 1 and 2 Respectively(m/s)

sigma_1 = 0.11;                           % Stage 1 structural mass fraction
sigma_2 = 0.11;                           % Stage 2 structural mass fraction
sigma_3 = 0.11;                           % Stage 3/0 structural mass fraction
isp_1 = 330;                              % Stage 1 Isp (s)
isp_2 = 380;                              % Stage 2 Isp (s)
isp_3 = 380;                              % Stage 3 Isp (s)

%==========================================================================
% code
for j=1:2                             % Mission Selection
        if j==1
            m_0_1 = [];
            DV_split = [];
            Zephyr_1DV = [];
            Zephyr_2DV = [];
            for DV_percent1 = 0.15:0.05:0.50
            for DV_percent3 = 0.15:0.05:0.50
                DV_percent2 = 1-(DV_percent1+DV_percent3);
                    
                V_stage3 = V_ideal(j+1)*DV_percent3;
                DV_req3 = V_stage3+V_Orbital(j+1);
                mu_3 = exp((DV_req3)/(isp_3*g0));
                m_p3 = m_payload(j+1)*(((mu_3-1)*(1-sigma_3))./(1-(sigma_3*mu_3)));
                m_s3 = m_p3*((sigma_3)/(1-sigma_3));
                m_03 = m_p3+m_s3+m_payload(j+1);
        
                V_stage2 = V_ideal(j+1)*DV_percent2;
                DV_req2 = V_stage2;
                mu_2 = exp((DV_req2)/(isp_2*g0));
                m_p2 = m_03.*(((mu_2-1)*(1-sigma_2))./(1-(sigma_2*mu_2)));
                m_s2 = m_p2*((sigma_2)/(1-sigma_2));
                m_02 = m_p2+m_s2;
            
                V_stage1 = V_ideal(j+1)*DV_percent1;
                DV_req1 = V_stage1+V_Recover+V_loss(j+1);
                mu_1 = exp((DV_req1)/(isp_1*g0));
                m_p1 = (m_02+m_03).*(((mu_1-1)*(1-sigma_1))./(1-(sigma_1*mu_1)));
                m_s1 = m_p1*((sigma_1)/(1-sigma_1));
                m_01 = m_p1+m_s1;
                    
                m_1i = m_01+m_02+m_payload(j);
                m_1f = m_s1+m_02+m_payload(j);
                DV_req1_1 = isp_1*g0*log(m_1i/m_1f);
                m_2i = m_02+m_payload(j);
                m_2f = m_s2+m_payload(j);
                DV_req2_1 = isp_2*g0*log(m_2i/m_2f);
                DV_req_1 = DV_req2_1+DV_req1_1;
                if DV_req_1 >= V_ideal(j)+V_loss(j)+V_Recover+V_Orbital(j)
                   m_0_1 = [m_0_1,m_1i];
                   DV_split = [DV_split,[DV_percent3;DV_percent2;DV_percent1]];
                end
                    
            end
            end
            [m_optimal,I] = min(m_0_1);
            Zephyr_mass = m_optimal;
            Zephyr_split = DV_split(:,I);
                
            DV_percent3 = Zephyr_split(1);
            DV_percent2 = Zephyr_split(2);
            DV_percent1 = Zephyr_split(3);
                
            V_stage3 = V_ideal(j+1)*DV_percent3;
            DV_req3 = V_stage3+V_Orbital(j+1);
            mu_3 = exp((DV_req3)/(isp_3*g0));
            m_p3 = m_payload(j+1)*(((mu_3-1)*(1-sigma_3))./(1-(sigma_3*mu_3)));
            m_s3 = m_p3*((sigma_3)/(1-sigma_3));
            m_03 = m_p3+m_s3+m_payload(j+1);
        
            V_stage2 = V_ideal(j+1)*DV_percent2;
            DV_req2 = V_stage2;
            mu_2 = exp((DV_req2)/(isp_2*g0));
            m_p2 = m_03.*(((mu_2-1)*(1-sigma_2))./(1-(sigma_2*mu_2)));
            m_s2 = m_p2*((sigma_2)/(1-sigma_2));
            m_02 = m_p2+m_s2;
        
            V_stage1 = V_ideal(j+1)*DV_percent1;
            DV_req1 = V_stage1+V_Recover+V_loss(j+1);
            mu_1 = exp((DV_req1)/(isp_1*g0));
            m_p1 = (m_02+m_03).*(((mu_1-1)*(1-sigma_1))./(1-(sigma_1*mu_1)));
            m_s1 = m_p1*((sigma_1)/(1-sigma_1));
            m_01 = m_p1+m_s1;
                
            m_0 = m_01+m_02+m_03;
                
            Zephyr_2DV = {'Zephyr Mission 2','Stage 3','Stage 2','Stage 1';'Sigma',sigma_3,sigma_2,sigma_1;'Isp',isp_3,isp_2,isp_1;'%DV Split',DV_percent3,DV_percent2,DV_percent1;'Propellant Mass',m_p3,m_p2,m_p1;'Structural mass',m_s3,m_s2,m_s1;'Total Stage Mass',m_03,m_02,m_01;'Total Mass',' ',' ',m_0};
            Zephyr2 = [sigma_3,sigma_2,sigma_1;isp_3,isp_2,isp_1;DV_percent3,DV_percent2,DV_percent1;m_p3,m_p2,m_p1;m_s3,m_s2,m_s1;m_03,m_02,m_01;0,0,m_0];
            writecell(Zephyr_2DV,'LVMasses/Zephyr-2MassEstimate.csv')
                
            m_p3 = 0;
            m_s3 = 0;
            m_03 = 0;
            DV_percent3 = 0;
                
            mu_2 = (m_p2+(m_payload(j)*(1-sigma_2)))/((m_p2*sigma_2)+(m_payload(j)*(1-sigma_2)));
            DV_req2 = log(mu_2)*isp_2*g0;
            V_stage2 = DV_req2-V_Orbital(j);
            DV_percent2 = V_stage2/V_ideal(j);
            m_02 = m_02+m_payload(j);
                
            mu_1 = (m_p1+(m_02*(1-sigma_1)))/((m_p1*sigma_1)+(m_02*(1-sigma_1)));
            DV_req1 = log(mu_1)*isp_1*g0;
            V_stage1 = DV_req1-V_Recover-V_loss(j);
            DV_percent1 = V_stage1/V_ideal(j);
                
            m_0 = m_01+m_02+m_payload(j);
            
            Zephyr_1DV = {'Zephyr Mission 1','Stage 3','Stage 2','Stage 1';'Sigma',sigma_3,sigma_2,sigma_1;'Isp',isp_3,isp_2,isp_1;'%DV Split',DV_percent3,DV_percent2,DV_percent1;'Propellant Mass',m_p3,m_p2,m_p1;'Structural mass',m_s3,m_s2,m_s1;'Total Stage Mass',m_03,m_02,m_01;'Total Mass',' ',' ',m_0};
            Zephyr1 = [sigma_3,sigma_2,sigma_1;isp_3,isp_2,isp_1;DV_percent3,DV_percent2,DV_percent1;m_p3,m_p2,m_p1;m_s3,m_s2,m_s1;m_03,m_02,m_01;0,0,m_0];
            writecell(Zephyr_1DV,'LVMasses/Zephyr-1MassEstimate.csv')
        end
end