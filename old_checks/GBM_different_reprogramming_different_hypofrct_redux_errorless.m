%% Basic 
clc;
clear;

load('fit_result_data_GBM.mat')

%% Defining Variables
% Radiotherapy model
par = parameters;
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [2 1]; % dose fraction sizes in Gy
Frac = [30 60]; 
len_D = length(Doses);
line_lengths = [0.8 2.3];
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 2500; % simulation length after last fraction

% ODE model

p = .505; % self renewal probability 

h = 10^(5); z = 1; % control the strenth of inhibitory signal

weak_feedback_Q = false;
if weak_feedback_Q
    l = 10^(-7); % weak feedback
else
    l = 10^(3); % strong feedback 
end
n = 1; low = -6; high = -2; c_l = high-low+1;
C = [0 5.196*10.^(linspace(low,high,c_l))]; 
color = cell(1,length(C));
for jj=1:c_l
    color{jj+1} = [jj/c_l (c_l-jj)/c_l 1];%blue to red [jj/c_l (c_l-jj)/c_l 1]; cyan to magenta
end
color{1} = [0 0 0]
% dedifferentiation rate
r1 = log(2)/DT; % growth rate of CSC
r2 = log(2)/DT; % growth rate of DCC
d  = log(2)/DT; % death rate of DCC

% Initial Condition
total_start_frac = 0.0005/64;
ROI_radius = 1; % Radius of the simulation region of intrest
rho = 10^9; % density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*pi()*ROI_radius^3*rho;
total_start_cell_num = total_cell_num*total_start_frac;
for jj=1:len_D
    x =[]; y1 = []; y2 = [];
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
            % Defining Variables
            D = Doses(jj);
            frac_num = Frac(jj);
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
            options = odeset('Refine',1e0, 'maxstep', 1e0);
            % With treatment 
            [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[0, treat_days(1)],[sc_start tc_start], options);
            [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[0, treat_days(1)],[sc_start tc_start], options);
            % Without treatment 
            [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[0, acq_days_after_RT],[sc_start tc_start], options);

            figure(1)
            hold on
            plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
            plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', line_lengths(jj)) 
            title('\fontsize{14} GBM Interaction')  
            xlabel('\fontsize{12} Time (Days)')
            ylabel('\fontsize{12} Survivall cell number (log)')
            axis([0 acq_days_after_RT -inf inf])
            hold off 

    %         figure(10)
    %         hold on
    %         plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
    %         plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', line_lengths(jj)) 
    %         xlabel('Time (Days)')
    %         ylabel({'Survivall cell'; 'number (log)'})
    %         axis([0 acq_days_after_RT -inf inf])
    %         hold off


            x = [T(end,1)];
            y1 = [U(end,1)*total_cell_num];
            y2 = [U(end,2)*total_cell_num];
            y3 = [U(end,1)];
            y4 = [U(end,2)];
            S = sum(U,2);
            y5 = [U(end,1)./S(end)*100];
            y6 = [p./(1+l.*(U(end,2)).^n)];
            y7 = [r1./(1+h.*(U(end,2)).^z)];


            for i = 1:length(sim_resume_days)

                if sum(U(:,2) < 0) > 0
                    break
                end
            %%%%%%% with stem cell %%%%%%%%%
                u = U(end,1);
                v = U(end,2);
                u_new = u*exp(-a1*D-b1*D^2) + c*v*D; % apply RT 
                v_new = v*exp(-a2*D-b2*D^2) - c*v*D; % apply RT

                [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n),[sim_resume_days(i), treat_days(i+1)],[u_new v_new]);
                x = [x T(1,1) T(end,1)];
                y1 = [y1 U(1,1)*total_cell_num U(end,1)*total_cell_num];
                y2 = [y2 U(1,2)*total_cell_num U(end,2)*total_cell_num];
                y3 = [y3 U(1,1) U(end,1)];
                y4 = [y4 U(1,2) U(end,2)]; 
                SS = sum(U,2);
                y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
                y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
                y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];

                figure(1)
                hold on
                plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', line_lengths(jj)) 
                plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj))
                plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', line_lengths(jj))
                hold off

    %             figure(10)
    %             hold on
    %             plot(Ta(:,1),log(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
    %             plot(Ta(:,1),log(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', line_lengths(jj)) 
    %             plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj))
    %             plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', line_lengths(jj))
    %             hold off



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
            plot(T(:,1),U(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
            plot(T(:,1),U(:,2)*total_cell_num,'color', color{k},'LineWidth', line_lengths(jj)) 
            hold off 

    %         figure(10)
    %         hold on
    %         plot(T(:,1),log(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
    %         plot(T(:,1),log(U(:,2)*total_cell_num),'color', color{k},'LineWidth', line_lengths(jj)) 
    %         hold off 



        end

        %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
    end
end
fig1 = figure(1);
ax1 = fig1.CurrentAxes;
%set(ax1,'YScale','log')
%set(ax1,'ylim',[exp(3) exp(14)])
if weak_feedback_Q
    fdbk_type = 'Weak';
else
    fdbk_type = 'Strong';
end
title(['\fontsize{14} ',fdbk_type,' Fdbk, Dose:',num2str(Doses(2)),', Days:',num2str(Frac(2))])
prefix = ['C:\Users\jhvo9\Documents\vojh MCSB Repo\Projects\Yu\Varying Fractionation\Better Visuals\',lower(fdbk_type),'_',num2str(Doses(2)),'_',num2str(Frac(2))];
% savefig([prefix '.fig'])
% saveas(gcf,[prefix '.png'])
%fig2 = figure(10);
%[h_m h_i]=inset_fractionation_schedules(fig1,fig2);
%set(h_i,'xtick',100:20:160,'xlim',[99,160])

