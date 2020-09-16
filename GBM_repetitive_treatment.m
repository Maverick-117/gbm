%% Basic 
clc;
clear;

load('fit_result_data_GBM.mat')

%% Defining Variables
% Radiotherapy model
par = parameters;
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
% Doses = [2 5 10]; % dose fraction sizes in Gy
% Frac = [30 10 3]; 
Doses = [1 2 4]; % dose fraction sizes in Gy
Frac = [65 30 13]; 
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 2500; % simulation length after last fraction

% ODE model
color = {'k','r','b'};
p = .505; % self renewal probability 
x =[]; y1 = []; y2 = [];
h = 10^5; z = 1; % control the strenth of inhibitory signal
l = 10^(-7); n = 1;
c = 5.196*10^(-3); % dedifferentiation rate
r1 = 0.1777; % growth rate of CSC
r2 = 0.1777; % growth rate of DCC
d  = 0.1777/5; % death rate of DCC

% Initial Condition
total_start_frac = 0.0005/64;
ROI_radius = 1; % Radius of the simulation region of intrest
rho = 10^9; % density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*pi()*ROI_radius^3*rho;
total_start_cell_num = total_cell_num*total_start_frac;

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
    for k = 1:length(Doses)
        % Defining Variables
        D = Doses(k);
        frac_num = Frac(k);
        sc_start = total_start_frac*F;
        tc_start = total_start_frac-sc_start;
        % Defining treatment days and ODE simulation start time after each fraction
        weeks = floor(frac_num/5);
        total_days = frac_num + 2*weeks;
        acq_end = treat_start + total_days + acq_days_after_RT - 1;
        A = repmat([1 1 1 1 1 0 0], 1, weeks+1);
        A_new = A(1:total_days);
        treat_days = find(A_new==1)+treat_start-1; 
        treat_days = [treat_days 400+treat_days 800+treat_days 1200+treat_days];
        sim_resume_days = treat_days+10/(60*24);
        treat_days = [treat_days acq_end];
        % ODE simulation before first fraction of RT
        options = odeset('Refine',1, 'maxstep', 1);
         % With treatment 
        [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[0, treat_days(1)],[sc_start tc_start], options);
        [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[0, treat_days(1)],[sc_start tc_start], options);
        % Without treatment 
        [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[0, acq_days_after_RT],[sc_start tc_start], options);
         
        
        figure(1)
        hold on
        plot(Tb(:,1),log(Ub(:,1)*total_cell_num),'g','LineStyle','--','LineWidth', 2) % stem cell
        plot(Tb(:,1),log(Ub(:,2)*total_cell_num),'g','LineWidth', 2) 
%       title({'\fontsize{14} GBM compartment interaction without feedback'})
%       title({'\fontsize{14} GBM compartment interaction';'\fontsize{12} with feedback on m_u and m_v'})  
%       title({'\fontsize{14} GBM compartment interaction';'\fontsize{12}with feedback on p'})           
        title({'\fontsize{14} GBM compartment interaction';'\fontsize{12} with feedback on m_u,m_v and p'})
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} Survivall cell number (log)')
        axis([0 acq_days_after_RT -inf inf])
        hold off 
        
        figure(10)
        hold on
        plot(Tb(:,1),log(Ub(:,1)*total_cell_num),'g','LineStyle','--','LineWidth', 2) % stem cell
        plot(Tb(:,1),log(Ub(:,2)*total_cell_num),'g','LineWidth', 2) 
        axis([0 acq_days_after_RT -inf inf])
        xlabel('Time (Days)')
        ylabel({'Survivall cell'; 'number (log)'})
        hold off
     
        figure(2)
        hold on
        plot(Tb(:,1),Ub(:,1)./sum(Ub,2)*100,'g','LineWidth', 2)  
%       title({'\fontsize{14} GBM compartment interaction without feedback'})
%       title({'\fontsize{14} GBM compartment interaction';'\fontsize{12} with feedback on m_u and m_v'})  
%       title({'\fontsize{14} GBM compartment interaction';'\fontsize{12}with feedback on p'})           
        title({'\fontsize{14} GBM compartment interaction';'\fontsize{12} with feedback on m_u,m_v and p'}) 
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} CSCs percentage (%)')
        axis([0 acq_days_after_RT -inf inf])
        hold off
        
        figure(11)
        hold on
        plot(Tb(:,1), Ub(:,1)./sum(Ub,2)*100,'g','LineWidth', 2) 
        axis([0 acq_days_after_RT -inf inf])
        xlabel('Time (Days)')
        ylabel({'CSC'; 'percentage (%)'})
        hold off
        
        figure(4)
        hold on
        axis([0 acq_days_after_RT -inf inf])
        title('\fontsize{14} Self Renewal Probability')
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} p value')
        plot(Tb(:,1),p./(1+l.*(Ub(:,2)).^n),'g','LineWidth', 2)
        hold off
        
        figure(5)
        hold on
        axis([0 acq_days_after_RT 0.025 inf])
        title('\fontsize{14} Division Rates')
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} m_u = m_v value')
        plot(Tb(:,1),r1./(1+h.*(Ub(:,2)).^z),'g','LineWidth', 2)
        hold off
        
        
        x = [T(end,1)];
        y1 = [log(U(end,1)*total_cell_num)];
        y2 = [log(U(end,2)*total_cell_num)];
        y3 = [U(end,1)];
        y4 = [U(end,2)];
        S = sum(U,2);
        y5 = [U(end,1)./S(end)*100];
        y6 = [p./(1+l.*(U(end,2)).^n)];
        y7 = [r1./(1+h.*(U(end,2)).^z)];

        
        for i = 1:length(sim_resume_days)
        %%%%%%% with stem cell %%%%%%%%%
            u = U(end,1);
            v = U(end,2);
            u_new = u*exp(-a1*D-b1*D^2) + c*v*D; % apply RT 
            v_new = v*exp(-a2*D-b2*D^2) - c*v*D; % apply RT
            
            [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[sim_resume_days(i), treat_days(i+1)],[u_new v_new]);
            x = [x T(1,1) T(end,1)];
            y1 = [y1 log(U(1,1)*total_cell_num) log(U(end,1)*total_cell_num)];
            y2 = [y2 log(U(1,2)*total_cell_num) log(U(end,2)*total_cell_num)];
            y3 = [y3 U(1,1) U(end,1)];
            y4 = [y4 U(1,2) U(end,2)]; 
            SS = sum(U,2);
            y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
            y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
            y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];
            
            figure(1)
            hold on
            plot(Ta(:,1),log(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),log(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', 2)
            plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', 2)
            hold off
            
            figure(10)
            hold on
            plot(Ta(:,1),log(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),log(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', 2)
            plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', 2)
            hold off
            
            
            figure(2)
            hold on
            plot(Ta(:,1),Ua(:,1)./sum(Ua,2)*100,'color',color{k},'LineWidth', 2) %somehow needed or not needed
            plot(x(1:end-1),y5(1:end-1),'color', color{k},'LineWidth', 2)
            hold off
            
            figure(11)
            hold on
            plot(Ta(:,1),Ua(:,1)./sum(Ua,2)*100,'color',color{k},'LineWidth', 2) %somehow needed or not needed
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
        
        x = [];
        y1 = [];
        y2 = [];
        y3 = [];
        y4 = [];
        y5 = [];
        y6 = [];
        y7 = [];
        
        figure(1)
        hold on
        plot(T(:,1),log(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
        plot(T(:,1),log(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
        hold off 
        
        figure(10)
        hold on
        plot(T(:,1),log(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
        plot(T(:,1),log(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
        hold off 
        
        figure(2)
        hold on
        plot(T(:,1),U(:,1)./sum(U,2)*100,'color', color{k},'LineWidth', 2) 
        hold off
        
        figure(11)
        hold on
        plot(T(:,1),U(:,1)./sum(U,2)*100,'color', color{k},'LineWidth', 2) 
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
        
    clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
end
        
fig1 = figure(1);
fig2 = figure(10);
[h_m h_i]=inset_rep(fig1,fig2);
set(h_i,'xtick',100:20:160,'xlim',[99,160])

fig3 = figure(2);
fig4 = figure(11);
[h_m h_i]=inset(fig3,fig4);
set(h_i,'xtick',100:20:160,'xlim',[99,160])
