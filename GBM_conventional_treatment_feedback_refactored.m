%% Basic 
clc;
clear;
close all
load('fit_result_data_GBM.mat')
sig = 10^16 * sig; 
%% Defining Variables and Initialization
% Radiotherapy model
par = parameters;%[.052 .01 .071 .197 .203];%   
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [2]; % dose fraction sizes in Gy
Frac = [30]; 
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
chi = 10^4; mu_bar = C(1); % chi too high causes negative values
% negative values seems like a general issue here...


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
suffix = '_zero_case';
saveQ = false; 

logQ = false; srvQ = true; logFracQ = false;
srv_start = 0; %initializing survivin amount
cont_p_a = 10^4; cont_p_b = 10*cont_p_a; compt_mult = 100;
if srvQ
    %srvn_csc = 3.6; srvn_dcc = 0.05;
    % model is more sensitive to this than to cont_p
    zeta_mult1 = 1;
    zeta_mult2 = 1; %; 3.6, 0.5; 3.6, 5];
else
    %srvn_csc = 0; srvn_dcc = 0;
    zeta_mult1 = 0;
    zeta_mult2 = 0;
end
srvn_zeta = [3.6 * zeta_mult1, 0.05 * zeta_mult2];
srvn_csc = srvn_zeta(1); srvn_dcc = srvn_zeta(2);
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
        D = Doses(g);
        frac_num = Frac(g);
        c = C(k);
        sc_start = total_start_frac*F;
        tc_start = total_start_frac-sc_start;
        % Defining treatment days and ODE simulation start time after each fraction
        weeks = floor(frac_num/5);
        total_days = frac_num + 2*weeks;
        acq_end = treat_start + 150;%treat_start + total_days + acq_days_after_RT - 1;
        A = repmat([1 1 1 1 1 0 0], 1, weeks+1);
        A_new = A(1:total_days);
        treat_days = find(A_new==1)+treat_start-1;
        sim_resume_days = treat_days+10/(60*24); % ODE simulation resumes 10 minutes after RT
        treat_days = [treat_days acq_end];     
        % ODE simulation before first fraction of RT
        options = odeset('Refine',1, 'maxstep', 1);
        % With treatment 
        [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[0, treat_days(1)],[sc_start tc_start srv_start], options);
        [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[0, treat_days(1)],[sc_start tc_start srv_start], options);
        % Without treatment 
        [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[0, acq_end],[sc_start tc_start srv_start], options);
        UU = U(:,3) .* exp(-sig * (T)); 
        UUa = Ua(:,3) .* exp(-sig * (Ta)); 
        UUb = Ub(:,3) .* exp(-sig * (Tb)); % using an integrating factor to bypass stiffness
        U(:,3) = UU; Ua(:,3) = UUa; Ub(:,3) = UUb;
        %% pre-therapy growth dynamics and plotting
        
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
        S = U(:,1)+U(:,2);
        y5 = [U(end,1)./S(end)*100];
        y6 = [p./(1+l.*(U(end,2)).^n)];
        y7 = [r1./(1+h.*(U(end,2)).^z)];

        %% radiotherapy dynamics and plotting
        for i = 1:length(sim_resume_days)
        %%%%%%% with stem cell %%%%%%%%%
            LQ_param = {a1, b1, a2, b2, c, D};
            [u_new,v_new, s_new,SF_U, SF_V] = radiotherapy(U, LQ_param, surv_vec);
            [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[sim_resume_days(i), treat_days(i+1)],[u_new v_new s_new]);
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
            SS = U(:,1)+U(:,2);
            y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
            y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
            y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];
            xx = x(1:end-1);
            yy = y6(1:end-1);
            yyy = y7(1:end-1);
           
        end
    end
        
    %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
end

% containers combining all similar vectors
u_sc  = [Ua(:,1)'*total_cell_num y1(2:end-2) U(:,1)'*total_cell_num];
u_dc  = [Ua(:,2)'*total_cell_num y2(2:end-2) U(:,2)'*total_cell_num];
u_srv = [zeros(size(Ua(:,1)')) ys(2:end-2) U(:,3)'];
t_vec = [Ta' x(2:end-2) T'];

fg_pops = figure()
hold on
plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', 4,'DisplayName','no treatment (CSC)') % stem cell
plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', 4,'DisplayName', 'no treatment (DCC)') 
plot(t_vec,u_sc,'color', color{k},'LineStyle','--','LineWidth', 2,'DisplayName', 'treatment (CSC)') % stem cell
plot(t_vec,u_dc,'color', color{k},'LineWidth', 2,'DisplayName', 'treatment (DCC)') 
if logQ
    set(gca, 'YScale', 'log')
    ylabel('\fontsize{12} Survival cell number (log10)')
else
    ylabel('\fontsize{12} Survival cell number')
end

if weak_feedback_Q
    title({'\fontsize{12} with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
else
    title({'\fontsize{12} with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
end

xlabel('\fontsize{12} Time (Days)')
axis([0 acq_end 7e3 2e7])
legend()
hold off 

fg_total = figure()
hold on
plot(Tb(:,1),(Ub(:,1)+Ub(:,2))*total_cell_num,'g','LineWidth', 4,'DisplayName','no treatment total')
plot(t_vec,u_sc+u_dc,color{k},'LineStyle','-.','LineWidth', 2,'DisplayName', 'treatment total')

if weak_feedback_Q
    title({'\fontsize{12} Total size with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
else
    title({'\fontsize{12} Total size with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
end

xlabel('\fontsize{12} Time (Days)')
axis([0 acq_end 7e3 2e7])
if logQ
    set(gca, 'YScale', 'log')
    ylabel('\fontsize{12} Total survival cell number (log10)')
else

    ylabel('\fontsize{12} Total survival cell number')
end
legend()
hold off 

% plotting survivin
fg_srv = figure()
hold on
if logQ
    set(gca, 'YScale', 'log')
    ylabel({'Survivin Fraction (log10)'})
else
    ylabel({'Survivin Fraction'})
end
plot(t_vec,u_srv,'color', color{k},'LineStyle',':','LineWidth', 2)
xlabel('Time (Days)')
title({'Survivin over time'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
axis([0 acq_end 0 2.5e-4])
hold off

fg_alph = figure()
hold on
plot(t_vec,1 ./(1+compt_mult*cont_p_a*u_srv*1),'color', color{k+1},'LineStyle',':','LineWidth', 2,'DisplayName','\alpha_U')
plot(t_vec,1 ./(1+cont_p_a*u_srv*1),'color', color{k+2},'LineStyle',':','LineWidth', 2,'DisplayName','\alpha_V')
xlabel('Time (Days)')
ylabel('Feedback on \alpha_U and \alpha_V')
title({'Feedback on \alpha_U and \alpha_V over time'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
axis([0 acq_end 0 1])
legend()
hold off

fg_bet = figure()
hold on
plot(t_vec,1 ./(1+compt_mult*cont_p_b*u_srv*1),'color', color{k+1},'LineStyle',':','LineWidth', 2,'DisplayName','\beta_U')
plot(t_vec,1 ./(1+cont_p_b*u_srv*1),'color', color{k+2},'LineStyle',':','LineWidth', 2,'DisplayName','\beta_V')
xlabel('Time (Days)')
ylabel('Feedback on \beta_U and \beta_V')
title({'Feedback on \beta_U and \beta_V over time'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
axis([0 acq_end 0 1])
legend()
hold off

fg_mu = figure()
hold on
plot(t_vec,chi*u_srv ./(1+chi*u_srv*1),'color', color{k},'LineStyle',':','LineWidth', 2)
xlabel('Time (Days)')
ylabel('Positive feedback on \mu')
title({'Relative rate of de-differentiation over time'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
axis([0 acq_end 0 1])
hold off

fg_csc = figure()
hold on
plot(Tb(:,1),Ub(:,1)./(Ub(:,1)+Ub(:,2))*100,'g','LineWidth', 4,'DisplayName','no treatment') 
plot(t_vec,u_sc./(u_sc+u_dc)*100,'color', color{k},'LineWidth', 2,'DisplayName', 'treatment') 
if weak_feedback_Q
    title({'\fontsize{12} with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
else
    title({'\fontsize{12} with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
end
xlabel('\fontsize{12} Time (Days)')
ylabel('\fontsize{12} CSC percentage (%)')
axis([0 acq_end .1 200])
set(gca,'yscale','log')
legend()
hold off

fg_prob = figure()
hold on
axis([0 acq_end 0 1])
if weak_feedback_Q
    title({'\fontsize{14} Self Renewal Probability';'\fontsize{12} feedback on m_u, m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
else
    title({'\fontsize{14} Self Renewal Probability';'\fontsize{12} feedback on m_u, m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
end

%        
xlabel('\fontsize{12} Time (Days)')
ylabel('\fontsize{12} p value')
plot(Tb(:,1),p./(1+l.*(Ub(:,2)).^n),'g','LineWidth', 4,'DisplayName','no treatment')
plot(t_vec,p./(1+l.*(u_dc/total_cell_num).^n),'color', color{k},'LineWidth', 2,'DisplayName','treatment') 
legend()
hold off

fg_div = figure()
hold on
axis([0 acq_end 0 r1])
if weak_feedback_Q
    title({'\fontsize{14} Division Rate';'\fontsize{12} feedback on m_u, m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
else
    title({'\fontsize{14} Division Rate';'\fontsize{12} feedback on m_u, m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
end
xlabel('\fontsize{12} Time (Days)')
ylabel('\fontsize{12} m_u = m_v value')
plot(Tb(:,1),r1./(1+h.*(Ub(:,2)).^z),'g','LineWidth', 4,'DisplayName','no treatment')
plot(t_vec,r1./(1+h.*(u_dc/total_cell_num).^z),'color', color{k},'LineWidth', 2,'DisplayName','treatment') 
fg_div.CurrentAxes.YScale = 'log';
legend()
hold off


test_x = t_vec(length(Ta)+1:end); test_y = u_srv(length(Ta)+1:end);
sum(test_x .* test_y)
trapz(test_x,test_y)

low_lim_log_s = floor(log10(min(test_y(test_y ~= 0))));
high_lim_log_s = ceil(log10(max(test_y(test_y ~= 0))));
yticks = unique(10.^floor(log10(logspace(low_lim_log_s,high_lim_log_s,8))));
fg_srv_inset = fg_srv;
[fg_srv fg_srv_inset] = inset(fg_srv, fg_srv_inset);
set(fg_srv_inset,'yscale','log','ytick',yticks,'ylim',[10^low_lim_log_s 10^high_lim_log_s],'xlim',[treat_start acq_end])
title(fg_srv_inset,'')
low_lim_log_mu = floor(log10(min(chi * test_y(test_y ~= 0)./(1+chi*test_y(test_y ~= 0)))));
high_lim_log_mu = ceil(log10(max(chi * test_y(test_y ~= 0)./(1+chi*test_y(test_y ~= 0)))));
yticks = unique(10.^floor(log10(logspace(low_lim_log_mu,high_lim_log_mu,8))));
fg_mu_inset = fg_mu;
[fg_mu fg_mu_inset] = inset(fg_mu, fg_mu_inset);
set(fg_mu_inset,'yscale','log','ytick',yticks,'ylim',[10^low_lim_log_mu 10^high_lim_log_mu],'xlim',[treat_start acq_end])
title(fg_mu_inset,'')
%% final plot embroidery and saving
fg_pops_inset = fg_pops;
[fg_pops fg_pops_inset]=inset(fg_pops,fg_pops_inset);
if logFracQ
    set(fg_pops,'yscale','log');
end
set(fg_pops_inset,'xtick',linspace(100,acq_end,4),'xlim',[100,acq_end],'yscale','log')
title(fg_pops_inset, '')
set(gca,'yscale','log')
legend(fg_pops,'Location','northwest')

%cab 3 4 5 2 6 7 100 102 103 104 105 106
fg_drvtv = figure(); plot(T(1:end-1),diff(U(:,1))*total_cell_num./diff(T),'ro:','LineWidth',2);
xlim([treat_start acq_end]);
ylim([0 inf])
ylabel({['Stem Cell Derivative'] ['(Forward Finite Difference)']}); xlabel('Days'); 
title({['rate of change after treatment for'], ['mu\_bar = ' num2str(mu_bar) '; chi = ' num2str(chi)]});

%close 1 3 6
if weak_feedback_Q
    fdbk_type = 'Weak';
else
    fdbk_type = 'Strong';
end
% figure(); plot(T, 1 ./ (1+cont_p_a*U(:,3)))
% figure(); plot(T, 1 ./ (1+compt_mult*cont_p_b*U(:,3)))
prefix = ['C:\Users\jhvo9\Google Drive (vojh1@uci.edu)\a PhD Projects\GBM Modeling\Survivin Models\Model 2 - Dedifferentiation Term\figs_diagnose_sigma\'];
% paramFile = [prefix 'params.txt'];
% if ~isfile(paramFile)
%     fileID = fopen(paramFile,'wt');
%     fprintf(fileID,['Feedback type: ' num2str(fdbk_type) '\n'...
%         'gamma_alpha_V: ' num2str(cont_p_a) '\n'...
%         'gamma_beta_V: ' num2str(cont_p_b) '\n'...
%         'chi: ' num2str(chi) '\n\n'...
%         'Not fixed:\n* reprogramming rate order of magnitude\n* Radiotherapy Schedule\n* Survivin Decay Rate\n* Compartment Multiplier\n'...
%         'File format: pwr_dose_frac_sig_compt_mult.'])
%     fclose(fileID)
% end
filename_prefix = [num2str(pwr),'_',num2str(Doses(1)),'_',num2str(Frac(1)),...
    '_',num2str(sig),'_',num2str(compt_mult),'_',num2str(srvn_csc),'_',num2str(srvn_dcc)];
if saveQ
    savefig(fg_total,[prefix 'total\' filename_prefix '_total' suffix '.fig'])
    saveas(fg_total,[prefix 'total\' filename_prefix '_total' suffix '.png'])

    savefig(fg_prob,[prefix 'prob\' filename_prefix '_prob' suffix '.fig'])
    saveas(fg_prob,[prefix 'prob\' filename_prefix '_prob' suffix '.png'])

    savefig(fg_div,[prefix 'div\' filename_prefix '_div' suffix '.fig'])
    saveas(fg_div,[prefix 'div\' filename_prefix '_div' suffix '.png'])

    savefig(fg_csc,[prefix 'csc\' filename_prefix '_csc' suffix '.fig'])
    saveas(fg_csc,[prefix 'csc\' filename_prefix '_csc' suffix '.png'])
    
    savefig(fg_alph,[prefix 'fdbk\alph\' filename_prefix '_alph' suffix '.fig'])
    saveas(fg_alph,[prefix 'fdbk\alph\' filename_prefix '_alph' suffix '.png'])
    
    savefig(fg_bet,[prefix 'fdbk\beta\' filename_prefix '_beta' suffix '.fig'])
    saveas(fg_bet,[prefix 'fdbk\beta\' filename_prefix '_beta' suffix '.png'])
    
    fg_pops_astr=ancestor(fg_pops,'figure'); fg_pops_n = fg_pops_astr.Number;
    fg_mu_astr=ancestor(fg_mu,'figure'); fg_mu_n = fg_mu_astr.Number;
    fg_srv_astr=ancestor(fg_srv,'figure'); fg_srv_n = fg_srv_astr.Number;
    
    savefig(figure(fg_pops_n),[prefix 'pops\' filename_prefix '_pops' suffix '.fig'])
    saveas(figure(fg_pops_n),[prefix 'pops\' filename_prefix '_pops' suffix '.png'])
    
    savefig(figure(fg_mu_n),[prefix 'fdbk\mu\' filename_prefix '_mu' suffix '.fig'])
    saveas(figure(fg_mu_n),[prefix 'fdbk\mu\' filename_prefix '_mu' suffix '.png'])
    
    savefig(figure(fg_srv_n),[prefix 'srv\' filename_prefix '_srv' suffix '.fig'])
    saveas(figure(fg_srv_n),[prefix 'srv\' filename_prefix '_srv' suffix '.png'])
    
    savefig(fg_drvtv,[prefix 'drvtv\' filename_prefix '_drvtv' suffix '.fig'])
    saveas(fg_drvtv,[prefix 'drvtv\' filename_prefix '_drvtv' suffix '.png'])
end