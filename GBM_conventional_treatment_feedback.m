%% Basic 
clc;
clear;
load('fit_result_data_GBM.mat')
sig = 0.1;
%% Defining Variables and Initialization
% Radiotherapy model
par = parameters;%[.052 .01 .071 .197 .203];%   
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [6]; % dose fraction sizes in Gy
Frac = [10]; 
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 110; % simulation length after last fraction

% ODE model
color = {'r','k','b'};
p = .505; % self renewal probability 
x =[]; y1 = []; y2 = [];
h = 10^(5); z = 1; % control the strength of inhibitory signal
weak_feedback_Q = true;
if weak_feedback_Q
    l = 10^(-7); % weak feedback
    fdbk_type = 'Weak';
else
    l = 10^(3); % strong feedback 
    fdbk_type = 'Strong';
end
n = 1; 
pwr = 3;
C = [5.196*10^(-pwr)];
% dedifferentiation rate
r1 = 0.1777; % growth rate of CSC
r2 = 0.1777; % growth rate of DCC
d  = 0.1777; % death rate of DCC

% Initial Condition
total_start_frac = 0.2/64; % Victoria: 0.2/64; Nayeon: 0.0005/64
ROI_radius = 1; % Radius of the simulation region of intrest
rho = 10^9; % density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*pi()*ROI_radius^3*rho;
total_start_cell_num = total_cell_num*total_start_frac;
suffix = '_alt';
saveQ = false; logQ = false; srvQ = true;
srv_start = 0; %initializing survivin amount
cont_p_a = 10^4; cont_p_b = 10*cont_p_a; compt_mult = 100;
if srvQ
    srvn_csc = 3.6; srvn_dcc = 0.05;
else
    srvn_csc = 0; srvn_dcc = 0;
end
surv_vec = {cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc}; %assuming these control parameters are constant





%% Tumor Growth ODE and radiotherapy simulation
for g = 1:length(cell_lines)
    F = par(g,1);
    a1 = par(g,2);  
    b1 = par(g,3); 
    a2 = par(g,4); 
    b2 = par(g,5); 
    a = par_ab(g,1); 
    b = par_ab(g,2);
    fprintf(['Cell Line: ' cell_lines{g} '\r'])
    fprintf(['a = ' num2str(a) ' b = ' num2str(b) '\r']);
    for k = 1:length(C)
        %% looping over reprogramming values
         % Defining Variables
        D = Doses;
        frac_num = Frac;
        c = C(k);
        sc_start = total_start_frac*F;
        tc_start = total_start_frac-sc_start;
        % Defining treatment days and ODE simulation start time after each fraction
        weeks = floor(frac_num/5);
        total_days = frac_num + 2*weeks;
        acq_end = treat_start + total_days + acq_days_after_RT - 1;
        A = repmat([1 1 1 1 1 0 0], 1, weeks+1);
        A_new = A(1:total_days);
        treat_days = find(A_new==1)+treat_start-1;
        sim_resume_days = treat_days+10/(60*24); % ODE simulation resumes 10 minutes after RT
        treat_days = [treat_days acq_end];     
        % ODE simulation before first fraction of RT
        options = odeset('Refine',1, 'maxstep', 1);
        % With treatment 
        [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, treat_days(1)],[sc_start tc_start srv_start], options);
        [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, treat_days(1)],[sc_start tc_start srv_start], options);
        % Without treatment 
        [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, acq_end],[sc_start tc_start srv_start], options);
        UU = U(:,3) .* exp(-sig * (T)); 
        UUa = Ua(:,3) .* exp(-sig * (Ta)); 
        UUb = Ub(:,3) .* exp(-sig * (Tb)); % using an integrating factor to bypass stiffness
        U(:,3) = UU; Ua(:,3) = UUa; Ub(:,3) = UUb;
        %% pre-therapy growth dynamics and plotting
        figure(1)
        hold on
        if logQ
            plot(Tb(:,1),log10(Ub(:,1)*total_cell_num),'g','LineStyle','--','LineWidth', 4) % stem cell
            plot(Tb(:,1),log10(Ub(:,2)*total_cell_num),'g','LineWidth', 4) 
            ylabel('\fontsize{12} Survival cell number (log10)')
        else
            plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', 2) % stem cell
            plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', 4) 
            ylabel('\fontsize{12} Survival cell number')
        end
        if weak_feedback_Q
            title({'\fontsize{12} with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{12} with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
        
        xlabel('\fontsize{12} Time (Days)')
        axis([0 acq_end -inf inf])
        hold off 
        
        figure(100)
        hold on
        if logQ
            plot(Tb(:,1),log10((Ub(:,1)+Ub(:,2))*total_cell_num),'g','LineWidth', 4) % stem cell
            plot(T(:,1),log10((U(:,1)+U(:,2))*total_cell_num),color{k},'LineStyle','-.','LineWidth', 2)
            ylabel('\fontsize{12} Survival cell number (log10)')
        else
            plot(Tb(:,1),(Ub(:,1)+Ub(:,2))*total_cell_num,'g','LineWidth', 4) % stem cell
            plot(T(:,1),(U(:,1)+U(:,2))*total_cell_num,color{k},'LineStyle','-.','LineWidth', 2)
            ylabel('\fontsize{12} Survival cell number')
        end
        if weak_feedback_Q
            title({'\fontsize{12} Total size with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{12} Total size with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
        
        xlabel('\fontsize{12} Time (Days)')
        axis([0 acq_end -inf inf])
        hold off 
        
        figure(10)
        hold on
        if logQ
            plot(Tb(:,1),log10(Ub(:,1)*total_cell_num),'g','LineStyle','-.','LineWidth', 4) % stem cell
            plot(Tb(:,1),log10(Ub(:,2)*total_cell_num),'g','LineWidth', 4) 
            plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            %plot(Ta(:,1),log10(Ua(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel({'Survival cell'; 'number (log10)'})
        else
            plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', 4) % stem cell
            plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', 4) 
            plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            %plot(Ta(:,1),Ua(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel({'Survival cell'; 'number'})
        end
        xlabel('Time (Days)')
        axis([0 acq_end -inf inf])
        hold off
        
        % plotting survivin
        figure(101)
        hold on
        if logQ
            plot(Ta(:,1),log10(UUa*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel({'Survivin Fraction (log10)'})
        else
            plot(Ta(:,1),UUa*1,'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel({'Survivin Fraction'})
        end
        xlabel('Time (Days)')
        title('Survivin over time')
        axis([0 acq_end -inf inf])
        hold off
        
        figure(102)
        hold on
        plot(Ta(:,1),1 ./(1+cont_p_a*UUa*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        xlabel('Time (Days)')
        ylabel('Feedback on \alpha_V')
        title('Feedback on \alpha_V over time')
        axis([0 acq_end 0 1])
        hold off
        figure(103)
        hold on
        plot(Ta(:,1),1 ./(1+cont_p_b*UUa*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        xlabel('Time (Days)')
        ylabel('Feedback on \beta_V')
        title('Feedback on \beta_V over time')
        axis([0 acq_end 0 1])
        hold off
        figure(104)
        hold on
        plot(Ta(:,1),1 ./(1+compt_mult*cont_p_b*UUa*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        xlabel('Time (Days)')
        ylabel('Feedback on \beta_U')
        title('Feedback on \beta_U over time')
        axis([0 acq_end 0 1])
        hold off
        figure(105)
        hold on
        plot(Ta(:,1),1 ./(1+compt_mult*cont_p_a*UUa*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        xlabel('Time (Days)')
        ylabel('Feedback on \alpha_U')
        title('Feedback on \alpha_U over time')
        axis([0 acq_end 0 1])
        hold off
        
        
        
        figure(2)
        hold on
        plot(Tb(:,1),Ub(:,1)./(Ub(:,1)+Ub(:,2))*100,'g','LineWidth', 4) 
        if weak_feedback_Q
            title({'\fontsize{12} with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{12} with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
         
%        
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} CSC percentage (%)')
        axis([0 acq_end .1 200])
        hold off
        
        figure(11)
        hold on
        plot(Tb(:,1), Ub(:,1)./(Ub(:,1)+Ub(:,2))*100,'g','LineWidth', 4) 
        axis([0 acq_end .1 200])
        xlabel('Time (Days)')
        ylabel({'CSC'; 'percentage (%)'})
        hold off
        
        figure(4)
        hold on
        axis([0 acq_end -inf inf])
        if weak_feedback_Q
            title({'\fontsize{14} Self Renewal Probability';'\fontsize{12} feedback on m_u, m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{14} Self Renewal Probability';'\fontsize{12} feedback on m_u, m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
         
%        
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} p value')
        plot(Tb(:,1),p./(1+l.*(Ub(:,2)).^n),'g','LineWidth', 4)
        hold off
        
        figure(5)
        hold on
        axis([0 acq_end -inf inf])
        if weak_feedback_Q
            title({'\fontsize{14} Division Rate';'\fontsize{12} feedback on m_u, m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{14} Division Rate';'\fontsize{12} feedback on m_u, m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} m_u = m_v value')
        plot(Tb(:,1),r1./(1+h.*(Ub(:,2)).^z),'g','LineWidth', 4)
        hold off
          
        
        x = [T(end,1)];
        if logQ
            y1 = [log10(U(end,1)*total_cell_num)];
            y2 = [log10(U(end,2)*total_cell_num)];
            ys = [log10(UU(end)*1)];
        else
            y1 = [U(end,1)*total_cell_num];
            y2 = [U(end,2)*total_cell_num];
            ys = [UU(end)*1];
        end
        y3 = [U(end,1)];
        y4 = [U(end,2)];
        S = U(:,1)+U(:,2);
        y5 = [U(end,1)./S(end)*100];
        y6 = [p./(1+l.*(U(end,2)).^n)];
        y7 = [r1./(1+h.*(U(end,2)).^z)];

        %% radiotherapy dynamics and plotting
        for i = 1:length(sim_resume_days)
        %%%%%%% with stem cell %%%%%%%%%
            LQ_param = {a1, b1, a2, b2, c, D};
            [u_new,v_new, s_new,SF_U, SF_V] = radiotherapy(U, LQ_param, surv_vec);

            [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[sim_resume_days(i), treat_days(i+1)],[u_new v_new s_new]);
            %fprintf('old UU: %f, new UU: %f', UU, U(end,3) .* exp(-sig * (T-treat_days(i)) )
            UU = U(:,3) .* exp(-sig * (T-treat_days(i))); % using an integrating factor to bypass stiffness
            U(:,3) = UU;
            x = [x T(1,1) T(end,1)];
            if logQ
                y1 = [y1 log10(U(1,1)*total_cell_num) log10(U(end,1)*total_cell_num)];
                y2 = [y2 log10(U(1,2)*total_cell_num) log10(U(end,2)*total_cell_num)];
                ys = [ys log10(UU(1)*1) log10(UU(end)*1)];
            else
                y1 = [y1 U(1,1)*total_cell_num U(end,1)*total_cell_num];
                y2 = [y2 U(1,2)*total_cell_num U(end,2)*total_cell_num];
                ys = [ys UU(1)*1 UU(end)*1];
            end
            y3 = [y3 U(1,1) U(end,1)];
            y4 = [y4 U(1,2) U(end,2)]; 
            SS = U(:,1)+U(:,2);
            y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
            y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
            y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];
            
            figure(1)
            hold on
            if logQ
                plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
                plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
                %plot(Ta(:,1),log10(Ua(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
                %plot(T(1:end-1),log10(UU(1:end-1)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
                plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
                %plot(T(1:end-1),UU(1:end-1)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
                %plot(Ta(:,1),Ua(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', 2)
            plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', 2)
            %plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
            hold off
            
            figure(100)
            hold on
            
            if logQ
                plot(x(1:end-1),log10(10.^(y1(1:end-1))+10.^(y2(1:end-1))),'color', color{k},'LineStyle','-.','LineWidth', 2)
                %plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(x(1:end-1),y1(1:end-1)+y2(1:end-1),'color', color{k},'LineStyle','-.','LineWidth', 2)
                %plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            hold off
            figure(101)
            hold on
            plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
            hold off
            
            figure(102)
            hold on
            if logQ
                plot(x(1:end-1),1 ./(1+cont_p_a*(10.^ys(1:end-1))*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(x(1:end-1),1 ./(1+cont_p_a*ys(1:end-1)*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            hold off
            figure(103)
            hold on
            if logQ
                plot(x(1:end-1),1 ./(1+cont_p_b*(10.^ys(1:end-1))*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(x(1:end-1),1 ./(1+cont_p_b*ys(1:end-1)*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            hold off
            figure(104)
            hold on
            if logQ
                plot(x(1:end-1),1 ./(1+compt_mult*cont_p_b*(10.^ys(1:end-1))*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(x(1:end-1),1 ./(1+compt_mult*cont_p_b*ys(1:end-1)*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            hold off
            figure(105)
            hold on
            if logQ
                plot(x(1:end-1),1 ./(1+compt_mult*cont_p_a*(10.^ys(1:end-1))*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(x(1:end-1),1 ./(1+compt_mult*cont_p_a*ys(1:end-1)*1),'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            hold off
            figure(10)
            hold on
            if logQ
                plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
                plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
                plot(Ta(:,1),log10(Ua(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
            else
                plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
                plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
                plot(Ta(:,1),Ua(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
            end
            plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', 2)
            plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', 2)
            hold off
            
            
            figure(2)
            hold on
            plot(Ta(:,1),Ua(:,1)./(Ua(:,1)+Ua(:,2))*100,'color',color{k},'LineWidth', 2) %somehow needed or not needed
            plot(x(1:end-1),y5(1:end-1),'color', color{k},'LineWidth', 2)
            hold off
            
            figure(11)
            hold on
            plot(Ta(:,1),Ua(:,1)./(Ua(:,1)+Ua(:,2))*100,'color',color{k},'LineWidth', 2) %somehow needed or not needed
            plot(x(1:end-1),y5(1:end-1),'color', color{k},'LineWidth', 2)
            hold off
            
            figure(4)
            hold on
            plot(Ta(:,1),p./(1+l.*(Ua(:,2)).^n),'color', color{k},'LineWidth', 2) 
            plot(x(1:end-1),y6(1:end-1),'color',color{k},'LineWidth', 2)
            hold off
            
            figure(5)
            hold on
            plot(Ta(:,1),r1./(1+h.*(Ua(:,2)).^z),'color', color{k},'LineWidth', 2) 
            plot(x(1:end-1),y7(1:end-1),'color',color{k},'LineWidth', 2)
            hold off
            
            xx = x(1:end-1);
            yy = y6(1:end-1);
            yyy = y7(1:end-1);
           
        end
        %% post-therapy dynamics and plotting
%         x = [];
%         y1 = [];
%         y2 = [];
%         y3 = [];
%         y4 = [];
%         y5 = [];
%         y6 = [];
%         y7 = [];
        
        
        figure(1)
        hold on
        
        if logQ
            plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            %plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
        else
            plot(T(:,1),U(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),U(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            %plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        hold off 
        figure(100)
        hold on
        
        if logQ
            plot(T(:,1),log10((U(:,1)+U(:,2))*total_cell_num),'color', color{k},'LineStyle','-.','LineWidth', 2) % stem cell
            %plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
         else
            plot(T(:,1),(U(:,1)+U(:,2))*total_cell_num,'color', color{k},'LineStyle','-.','LineWidth', 2) % stem cell
            %plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        hold off 
        figure(101)
        hold on
        
        if logQ
            plot(T(:,1),log10(UU*1),'color', color{k},'LineStyle',':','LineWidth', 2)
         else
            plot(T(:,1),UU*1,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        hold off 
        
        figure(102)
        hold on
        plot(T(:,1),1./(1+cont_p_a*UU*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        hold off
        figure(103)
        hold on
        plot(T(:,1),1./(1+cont_p_b*UU*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        hold off
        figure(104)
        hold on
        plot(T(:,1),1./(1+compt_mult*cont_p_b*UU*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        hold off
        figure(105)
        hold on
        plot(T(:,1),1./(1+compt_mult*cont_p_a*UU*1),'color', color{k},'LineStyle',':','LineWidth', 2)
        hold off
        
        figure(10)
        hold on
        if logQ
            plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            %plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
        else
            plot(T(:,1),U(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),U(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            %plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        hold off 
        
        figure(2)
        hold on
        plot(T(:,1),U(:,1)./(U(:,1)+U(:,2))*100,'color', color{k},'LineWidth', 2) 
        hold off
        
        figure(11)
        hold on
        plot(T(:,1),U(:,1)./(U(:,1)+U(:,2))*100,'color', color{k},'LineWidth', 2) 
        hold off
        
        figure(4)
        hold on
        tt = T(:,1);
        pp = p./(1+l.*(U(:,2)).^n);
        plot([xx(end),tt(1)],[yy(end),pp(1)],'color',color{k},'LineWidth', 2)
        plot(T(:,1),p./(1+l.*(U(:,2)).^n),'color',color{k},'LineWidth', 2)
        hold off
       
        
        figure(5)
        hold on
        ttt = T(:,1);
        rr = r1./(1+h.*(U(:,2)).^z);
        plot([xx(end),ttt(1)],[yyy(end),rr(1)],'color',color{k},'LineWidth', 2)
        plot(T(:,1),r1./(1+h.*(U(:,2)).^z),'color',color{k},'LineWidth', 2)
        hold off
      
        
    end
        
    %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
end
%% final plot embroidery and saving
fig1 = figure(1);
fig2 = figure(10);
[h_m h_i]=inset(fig1,fig2);
set(h_i,'xtick',100:20:160,'xlim',[99,160])

fig3 = figure(2);
fig4 = figure(11);
[h_m h_i]=inset(fig3,fig4);
set(h_i,'xtick',100:20:160,'xlim',[99,160],'yscale','log')
set(h_m,'yscale','log')

cab 3 4 5 6 100 101 102 103 104 105
fig1 = figure(3);
fig2 = figure(4);
fig3 = figure(5);
fig4 = figure(6);

if weak_feedback_Q
    fdbk_type = 'Weak';
else
    fdbk_type = 'Strong';
end
% figure(); plot(T, 1 ./ (1+cont_p_a*U(:,3)))
% figure(); plot(T, 1 ./ (1+compt_mult*cont_p_b*U(:,3)))
prefix = ['C:\Users\jhvo9\Documents\vojh MCSB Repo\Projects\Yu\Varying Fractionation\Threshold\'];
filename_prefix = [lower(fdbk_type),num2str(pwr),'_',num2str(Doses(1)),'_',num2str(Frac(1))];
if saveQ
    savefig(fig1,[prefix 'total\' filename_prefix '_total' suffix '.fig'])
    saveas(fig1,[prefix 'total\' filename_prefix '_total' suffix '.png'])

    savefig(fig2,[prefix 'prob\' filename_prefix '_prob' suffix '.fig'])
    saveas(fig2,[prefix 'prob\' filename_prefix '_prob' suffix '.png'])

    savefig(fig3,[prefix 'div\' filename_prefix '_div' suffix '.fig'])
    saveas(fig3,[prefix 'div\' filename_prefix '_div' suffix '.png'])

    savefig(fig4,[prefix 'csc\' filename_prefix '_csc' suffix '.fig'])
    saveas(fig4,[prefix 'csc\' filename_prefix '_csc' suffix '.png'])
end

test_x = [x(1:end-2) T']; test_y = [ys(1:end-2) U(:,3)'];
sum(test_x .* test_y)
trapz(test_x,test_y)