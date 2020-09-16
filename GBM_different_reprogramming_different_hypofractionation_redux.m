%% Basic 
clc;
clear;

load('fit_result_data_GBM.mat')
warning_counter = 0; culprit1 = []; culprit2 = []; d_type = []; c_type = []; xc = [];
%% Defining Variables
% Radiotherapy model
par = parameters;
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [2 15]; % dose fraction sizes in Gy
Frac = [30 4]; 
len_D = length(Doses);
line_lengths = [0.8 2.3];
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 600; % simulation length after last fraction

% ODE model

p = .505; % self renewal probability 

h = 10^(5); z = 1; % control the strenth of inhibitory signal


n = 1; low = -5; high = -2; c_l = high-low+1;
C = [0 5.196*10.^(linspace(low,high,c_l))]; 
color = cell(1,length(C));
saveQ = true; logQ = true; srvQ = true;
for jj=1:c_l
    color{jj+1} = [jj/c_l (c_l-jj)/c_l 1];
end
color{1} = [0 0 0];
% dedifferentiation rate
r1 = log10(2)/DT; % growth rate of CSC
r2 = log10(2)/DT; % growth ra te of DCC
d  = log10(2)/DT; % death rate of DCC

% Initial Condition
total_start_frac = 0.0005/64;
ROI_radius = 1; % Radius of the simulation region of intrest
rho = 10^9; % density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*pi()*ROI_radius^3*rho;
total_start_cell_num = total_cell_num*total_start_frac;
ll = [false true];
suffix = '_alt';
cont_p_a = 1; cont_p_b = 0.1; compt_mult = 10;
if srvQ
    srvn_csc = 3.6; srvn_dcc = 0.05;
else
    srvn_csc = 0; srvn_dcc = 0;
end
surv_vec = {cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc}; %assuming these control parameters are constant

for ii=1:length(ll)
    weak_feedback_Q = ll(ii);
    if weak_feedback_Q
        l = 10^(-7); % weak feedback
        fdbk_type = 'Weak';
    else
        l = 10^(3); % strong feedback 
        fdbk_type = 'Strong';
    end
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
                srv_start = 0; %initializing survivin amount
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
                options = odeset('Refine',1, 'maxstep', 1);
                % With treatment 
                [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, treat_days(1)],[sc_start tc_start srv_start], options);
                [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, treat_days(1)],[sc_start tc_start srv_start], options);
                % Without treatment 
                [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, acq_days_after_RT],[sc_start tc_start srv_start], options);
                figure(ii)
                hold on

                if logQ
                    plot(Tb(:,1),log10(Ub(:,1)*total_cell_num),'g','LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                    plot(Tb(:,1),log10(Ub(:,2)*total_cell_num),'g','LineWidth', line_lengths(jj))
                    ylabel('\fontsize{12} Survival cell number (log)')
                else
                    plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                    plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', line_lengths(jj)) 
                    ylabel('\fontsize{12} Survival cell number')
                end
                title('\fontsize{14} GBM Interaction')  
                xlabel('\fontsize{12} Time (Days)')
                axis([0 acq_days_after_RT -inf inf])
                hold off 



                x = [T(end,1)];
                if logQ
                    y1 = [log10(U(end,1)*total_cell_num)];
                    y2 = [log10(U(end,2)*total_cell_num)];
                    ys = [log10(U(end,3)*total_cell_num)];
                else
                    y1 = [U(end,1)*total_cell_num];
                    y2 = [U(end,2)*total_cell_num];
                    ys = [U(end,3)*total_cell_num];
                end
                y3 = [U(end,1)];
                y4 = [U(end,2)];
                S = sum(U,2);
                y5 = [U(end,1)./S(end)*100];
                y6 = [p./(1+l.*(U(end,2)).^n)];
                y7 = [r1./(1+h.*(U(end,2)).^z)];


                for i = 1:length(sim_resume_days)
                %%%%%%% with stem cell %%%%%%%%%
                    LQ_param = {a1, b1, a2, b2, c, D};
                    [u_new,v_new, s_new] = radiotherapy(U,LQ_param, surv_vec);

                    [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[sim_resume_days(i), treat_days(i+1)],[u_new v_new s_new]);
                    x = [x T(1,1) T(end,1)];
                    if logQ
                        y1 = [y1 log10(U(1,1)*total_cell_num) log10(U(end,1)*total_cell_num)];
                        y2 = [y2 log10(U(1,2)*total_cell_num) log10(U(end,2)*total_cell_num)];
                        ys = [ys log10(U(1,3)*total_cell_num) log10(U(end,3)*total_cell_num)];
                    else
                        y1 = [y1 U(1,1)*total_cell_num U(end,1)*total_cell_num];
                        y2 = [y2 U(1,2)*total_cell_num U(end,2)*total_cell_num];
                        ys = [ys U(1,3)*total_cell_num U(end,3)*total_cell_num];
                    end
                    y3 = [y3 U(1,1) U(end,1)];
                    y4 = [y4 U(1,2) U(end,2)]; 
                    SS = sum(U,2);
                    y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
                    y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
                    y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];

                    figure(ii)
                    hold on
                    if logQ
                        plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                        plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
                        %plot(T(1:end-1),log10(U(1:end-1,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
                    else
                        plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                        plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
                        %plot(T(1:end-1),U(1:end-1,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
                    end
                    plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj))
                    lastwarn('', '');
                    plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', line_lengths(jj))
                    [warnMsg, warnId] = lastwarn();
                    if strcmp(warnId,'MATLAB:plot:IgnoreImaginaryXYPart')
                        fprintf("bananas.")
                        xc = [xc T(1,1) T(end,1)];
                        culprit2 = [culprit2 U(1,2) U(end,2)];
                        d_type = [d_type jj];
                        c_type = [c_type k];
                        warning_counter = warning_counter + 1;
                    end
                    hold off

        %             figure(10)
        %             hold on
        %             plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
        %             plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', line_lengths(jj)) 
        %             plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj))
        %             plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', line_lengths(jj))
        %             hold off



                    xx = x(1:end-1);
                    yy = y6(1:end-1);
                    yyy = y7(1:end-1);



                end

    %            x = [];
    %             y1 = [];
    %             y2 = [];
                y3 = [];
                y4 = [];
                y5 = [];
                y6 = [];
                y7 = [];

                figure(ii)
                hold on
                if logQ
                    plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                    plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', line_lengths(jj)) 
                    %plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
                else
                    plot(T(:,1),U(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
                    plot(T(:,1),U(:,2)*total_cell_num,'color', color{k},'LineWidth', line_lengths(jj)) 
                    %plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
                end
                hold off 

        %         figure(10)
        %         hold on
        %         plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', line_lengths(jj)) % stem cell
        %         plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', line_lengths(jj)) 
        %         hold off 



            end

            %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
        end
    end
    fig1 = figure(ii);
    ax1 = fig1.CurrentAxes;
    if logQ
        set(ax1,'ylim',[-3 6])
    else
        set(ax1,'ylim',[10^(-3) 10^(6)])
    end
    if weak_feedback_Q
        fdbk_type = 'Weak';
    else
        fdbk_type = 'Strong';
    end
    title(['\fontsize{14} ',fdbk_type,' Fdbk, Dose:',num2str(Doses(2)),' Gy, Days:',num2str(Frac(2)) ', Srvn Params: (' num2str(cont_p_a) ',' num2str(cont_p_b) ',' num2str(compt_mult) ')'])
    prefix = ['C:\Users\jhvo9\Documents\vojh MCSB Repo\Projects\Yu\Varying Fractionation\Better Visuals\',lower(fdbk_type),'_',num2str(Doses(2)),'_',num2str(Frac(2))];
    if saveQ
        savefig([prefix suffix '.fig'])
        saveas(gcf,[prefix suffix '.png'])
    end
    %close(gcf)
end

%fig2 = figure(10);
%[h_m h_i]=inset_fractionation_schedules(fig1,fig2);
%set(h_i,'xtick',100:20:160,'xlim',[99,160])

