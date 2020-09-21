%% Basic 
clc;
clear;

load('fit_result_data_GBM.mat')
warning_counter = 0; culprit1 = []; culprit2 = []; d_type = []; c_type = []; xc = [];
%% Defining Variables
% Radiotherapy model
par = parameters;
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % dose fraction sizes in Gy
Frac = [60 30 20 15 12 10 8 7 6 6 5 5 4 4 4 3 3 3 3 3]; 
len_D = length(Doses);
%sample_times = [1 5 10 20 30]; 
pre_sample_times = [100 300]; %600 900 1200 1500];
len_samples = length(pre_sample_times);

dose1offset = 0;
line_lengths = [0.8 2.3];
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 2500; % simulation length after last fraction
% ODE model

p = .505; % self renewal probability 

h = 10^(5); z = 1; % control the strenth of inhibitory signal


saveQ = false; logQ = true; weak_feedback_Q = true; srvQ = true;



if weak_feedback_Q
    l = 10^(-7); % weak feedback
else
    l = 10^(3); % strong feedback 
end
n = 1; 
%low = -3; high = -3; c_l = 1*(high-low+1);
%C = [0 5.196*10.^(linspace(low,high,c_l))]; % ramp up by factors of 2 later
C = [5.196*10^(-3)];
% [0 5.196*2.^(linspace(low,high,c_l))];
len_C = length(C);
times = zeros(len_D,len_samples);



% dedifferentiation rate
r1 = log10(2)/DT; % growth rate of CSC
r2 = log10(2)/DT; % growth rate of DCC
d  = log10(2)/DT; % death rate of DCC

% Initial Condition
total_start_frac = 0.2/64;
ROI_radius = 1; % Radius of the simulation region of intrest
rho = 10^9; % density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*pi()*ROI_radius^3*rho;
total_start_cell_num = total_cell_num*total_start_frac;

if srvQ
    srvn_csc = 3.6; srvn_dcc = 0.05;
    srvn_zeta = [0, 0; 1, 0.01;1, 0.1; 1, 1; 1, 10; 1, 100; 1, 1000];%[0, 0; 0.1, 0.01; 1, 0.1; 10, 1000; 1, 100];
    %[0, 0; 1, 0.01;1, 0.1; 1, 1; 1, 10; 1, 100; 1, 1000];
    srvn_zeta(:,1) = srvn_zeta(:,1) * srvn_csc;
    srvn_zeta(:,2) = srvn_zeta(:,2) * srvn_dcc;
    % model is more sensitive to this than to cont_p
else
    %srvn_csc = 0; srvn_dcc = 0;
    srvn_zeta = [0, 0];
end


cont_p_a = 10^4;
%cont_p = [0 0; 1 0.1; 1 1; 1 10; 1 100] * cont_p_a;
% the simulation seems to be sort of sensitive to control parameters
cont_p = [1 10] * cont_p_a;
%cont_p_a = 10^4; % \gamma_{\alpha_V}
%cont_p_b = [0.1 10]; % \gamma_{\beta_V}
compt_mult = 2; % \zeta_U; in reality, this should be higher than 1 since stem cells 
                % tend to more difficult to kill, so a higher denominator is more attractive
len_surv_cont = size(cont_p); 
len_surv_cont = len_surv_cont(1);
len_zeta = size(srvn_zeta); len_zeta = len_zeta(1);
dose_reprog_time_cont = zeros(len_D,len_C,len_samples,len_zeta,len_surv_cont);
N_c = len_zeta;
color = cell(1,N_c);
for jj=2:N_c
    color{jj} = [jj/N_c (N_c-jj)/N_c 1];
end
color{1} = [0 0 0];

for ii=1:len_surv_cont
    for jj=1:len_D
        x =[]; y1 = []; y2 = [];    %% Tumor Growth ODE and radiotherapy simulation
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
            for gg=1:len_zeta
                surv_vec = {cont_p(ii,1), cont_p(ii,2), compt_mult, srvn_zeta(gg,1), srvn_zeta(gg,2)}; %assuming these control parameters are constantd
                for k = 1:len_C
                    % Defining Variables
                    srv_start = 0; %initializing survivin amount
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
                    [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[0, acq_end],[sc_start tc_start srv_start], options);
                    U(:,3) = U(:,3) .* exp(-sig * (T)); 
                    Ua(:,3) = Ua(:,3) .* exp(-sig * (Ta)); 
                    Ub(:,3) = Ub(:,3) .* exp(-sig * (Tb)); % using an integrating factor to bypass stiffness


                    x = [T(end,1)];
                    if logQ
                        y1 = [log10(U(end,1)*total_cell_num)];
                        y2 = [log10(U(end,2)*total_cell_num)];
                        ys = [log10(U(end,3)*1)];
                    else
                        y1 = [U(end,1)*total_cell_num];
                        y2 = [U(end,2)*total_cell_num];
                        ys = [U(end,3)*1];
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
                        U(:,3) = U(:,3) .* exp(-sig * (T-treat_days(i)));
                        x = [x T(1,1) T(end,1)];
                        if logQ
                            y1 = [y1 log10(U(1,1)*total_cell_num) log10(U(end,1)*total_cell_num)];
                            y2 = [y2 log10(U(1,2)*total_cell_num) log10(U(end,2)*total_cell_num)];
                            ys = [ys log10(U(1,3)*1) log10(U(end,3)*1)];
                        else
                            y1 = [y1 U(1,1)*total_cell_num U(end,1)*total_cell_num];
                            y2 = [y2 U(1,2)*total_cell_num U(end,2)*total_cell_num];
                            ys = [ys U(1,3)*1 U(end,3)*1];
                        end
                        y3 = [y3 U(1,1) U(end,1)];
                        y4 = [y4 U(1,2) U(end,2)]; 
                        SS = sum(U,2);
                        y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
                        y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
                        y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];




                        xx = x(1:end-1);
                        yy = y6(1:end-1);
                        yyy = y7(1:end-1);



                    end

                    sample_times = 1:len_samples;
                    for i=1:len_samples
                        sample_times(i) = find(T < treat_start + pre_sample_times(i),1,'last');
                    end
                    if jj > 1
                        dose_reprog_time_cont(jj,k,:,gg,ii) = U(sample_times,1);
                    else
                        dose_reprog_time_cont(jj,k,:,gg,ii) = U(sample_times+dose1offset,1);
                    end
        %            x = [];
        %             y1 = [];
        %             y2 = [];
                    y3 = [];
                    y4 = [];
                    y5 = [];
                    y6 = [];
                    y7 = [];
                    if jj > 1
                        times(jj,:) = T(sample_times);
                    else
                        times(jj,:) = T([sample_times(1)+dose1offset sample_times(2:end)]); %same as before but with offset at first time to let times line up across dosage.
                    end

                
                end
            end
            %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
        end
    end
end
str_leg = arrayfun(@(x) num2str(x),C,'UniformOutput',false);
if weak_feedback_Q
    fdbk_type = 'Weak';
else
    fdbk_type = 'Strong';
end
s_ylabel = '#CSC/#CSC(conventional)';
s_xlabel = 'Dose (Gy)';
s_legend_title = 'Survivin Release Coefficients:';%'Survivin Control Params:';
save_prefix = ['C:\Users\jhvo9\Documents\vojh MCSB Repo\Projects\Yu\Varying Fractionation\Better Visuals\',lower(fdbk_type), '_relCSC_time_'];
xlim_vec = [min(Doses) max(Doses)];
ylim_vec = [0 1.35];
leg_pos = 'north';
line_types = cell(1,len_surv_cont);
for k=2:len_surv_cont
    line_types{k} = '--';
end
line_types{1} = '-';
suffix = '_alt'; % i.e. '_testN'
per_ratio = zeros(len_D,len_C,len_samples,len_surv_cont);
fig_start = 50;
for ii=1:len_surv_cont
    title_prefix = [fdbk_type,' Fdbk;(' num2str(cont_p(ii,1)) ',' num2str(cont_p(ii,2)) ',' num2str(compt_mult) '); Relative CSC at time ~ '];%[fdbk_type,' Fdbk; Relative CSC at time ~ '];
    for t=1:len_samples
        figX = figure(fig_start*ii+(t-1));
        hold on;
        for k=1:len_C
            for gg=1:len_zeta
                plot(Doses,dose_reprog_time_cont(:,k,t,gg,ii)/dose_reprog_time_cont(2,k,t,gg,ii),'o',...
                    'color', color{gg},'LineStyle',line_types{ii},'LineWidth',line_lengths(2),...
                    'DisplayName',['(' num2str(srvn_zeta(gg,1)) ',' num2str(srvn_zeta(gg,2)) ')']);
                % ['(' num2str(cont_p(ii,1)) ',' num2str(cont_p(ii,2)) ',' num2str(compt_mult) ')']
                per_ratio(:,k,t,gg,ii) = max(0,1-dose_reprog_time_cont(:,k,t,gg,ii)/dose_reprog_time_cont(2,k,t,gg,ii));
            end
        end
        lgdX = legend('Location',leg_pos);
        title(lgdX,s_legend_title)
        title([title_prefix,num2str(round(mean(times(:,t)),ceil(log10(mean(times(:,t)))),'significant')),' days'])
        ylabel(s_ylabel)
        xlabel(s_xlabel)
        xlim(xlim_vec)
        ylim(ylim_vec)
        hold off;
        if saveQ
            savefig([save_prefix num2str(sample_times(t)) suffix '.fig'])
            saveas(gcf,[save_prefix num2str(sample_times(t)) suffix '.png'])
        end
    end
end



per_ratio(per_ratio == 0) = nan;
                    

% fig_start2 = 60;
% control_ratio = cont_p(ii,2)/cont_p(ii,1);
% for t=1:len_samples
%     figX = figure(fig_start2*k+t-1);
%     hold on;
%     for k=1:len_C
%         for ii=1:len_surv_cont
%             scatter(sort(repmat(control_ratio(ii),1,length(Doses))),Doses,...
%                 100*per_ratio(:, k,t,ii),color{k})
%             title(['time: ~' num2str(mean(times(:,t))) ]);
%             xlabel('survivin control ratio');
%             ylabel('Doses')
%         end
%     end
%     hold off;
% end
% fig2 = figure(30);
% hold on;
% for k=1:len_C
%     plot(Doses,dose_reprog_time(:,k,1)/dose_reprog_time(2,k,1),'o',...
%         'color', color{k},'LineStyle',line_types{k},'LineWidth',line_lengths(2));
% end
% lgd2 = legend(strC,'Location',leg_pos);
% title(lgd2,s_legend_title)
% title([title_prefix,num2str(round(mean(times(:,1)),ceil(log10(mean(times(:,1)))),'significant')),' days'])
% ylabel(s_ylabel)
% xlabel(s_xlabel)
% xlim(xlim_vec)
% ylim(ylim_vec)
% hold off;
% if saveQ
%     savefig([save_prefix '1_' suffix '.fig'])
%     saveas(gcf,[save_prefix '1_' suffix '.png'])
% end

% % The figures below produce the same stuff for weak feedback
% fig3 = figure(31);
% hold on;
% for k=1:len_C
%     plot(Doses,dose_reprog_time(:,k,2)/dose_reprog_time(2,k,2),'o',...
%         'color', color{k},'LineStyle',line_types{k},'LineWidth',line_lengths(2));
% end
% lgd3 = legend(strC,'Location',leg_pos);
% title(lgd3,s_legend_title)
% title([title_prefix,num2str(round(mean(times(:,2)),ceil(log10(mean(times(:,2)))),'significant')),' days'])
% ylabel(s_ylabel)
% xlabel(s_xlabel)
% xlim(xlim_vec)
% ylim(ylim_vec)
% hold off;
% if saveQ
%     savefig([save_prefix '2_' suffix '.fig'])
%     saveas(gcf,[save_prefix '2_' suffix '.png'])
% end
% fig4 = figure(32);
% hold on;
% for k=1:len_C
%     plot(Doses,dose_reprog_time(:,k,3)/dose_reprog_time(2,k,3),'o',...
%         'color', color{k},'LineStyle',line_types{k},'LineWidth',line_lengths(2));
% end
% lgd4 = legend(strC,'Location',leg_pos);
% title(lgd4,s_legend_title)
% title([title_prefix,num2str(round(mean(times(:,3)),ceil(log10(mean(times(:,3)))),'significant')),' days'])
% ylabel(s_ylabel)
% xlabel(s_xlabel)
% xlim(xlim_vec)
% ylim(ylim_vec)
% hold off;
% if saveQ
%     savefig([save_prefix '3_' suffix '.fig'])
%     saveas(gcf,[save_prefix '3_' suffix '.png'])
% end
