%% Basic 
clc;
clear;

load('../fit_result_data_GBM.mat')
warning_counter = 0; culprit1 = []; culprit2 = []; d_type = []; c_type = []; xc = [];
%% Defining Variables
% Radiotherapy model
par = parameters;
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 60]; % dose fraction sizes in Gy
Frac = [60 30 20 15 12 10 8 7 6 6 5 5 4 4 4 3 3 3 3 3 1]; 
len_D = length(Doses);
%sample_times = [1 5 10 20 30]; len_samples = length(sample_times);
pre_sample_times = [100 300 600 900 1200 1500];
len_samples = length(pre_sample_times);
dose1offset = 0;
line_lengths = [0.8 2.3];
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 2500; % simulation length after last fraction
% ODE model

p = .505; % self renewal probability 

h = 10^(5); z = 1; % control the strenth of inhibitory signal

saveQ = false; logQ = true; weak_feedback_Q = true;
if weak_feedback_Q
    l = 10^(-7); % weak feedback
else
    l = 10^(3); % strong feedback 
end
n = 1; low = -3; high = -3; c_l = 1*(high-low+1);
C = [0 5.196*10.^(linspace(low,high,c_l))]; % ramp up by factors of 2 later
% [0 5.196*2.^(linspace(low,high,c_l))];
len_C = length(C);
dose_reprog_time = zeros(len_D,len_C,len_samples);
times = zeros(len_D,len_samples);

color = cell(1,len_C);
for jj=1:c_l
    color{jj+1} = [jj/c_l (c_l-jj)/c_l 1];
end
color{1} = [0 0 0]

% dedifferentiation rate
r1 = log(2)/DT; % growth rate of CSC
r2 = log(2)/DT; % growth ra te of DCC
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
            options = odeset('Refine',1, 'maxstep', 1);
            % With treatment 
            [T,U] = ode45(@(t,U) stem_ODE_feedback_no_srvn(t, U, r1, r2, d, p, h, z, l, n),[0, treat_days(1)],[sc_start tc_start], options);
            [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback_no_srvn(t, U, r1, r2, d, p, h, z, l, n),[0, treat_days(1)],[sc_start tc_start], options);
            % Without treatment 
            [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback_no_srvn(t, U, r1, r2, d, p, h, z, l, n),[0, acq_days_after_RT],[sc_start tc_start], options);



            x = [T(end,1)];
            if logQ
                y1 = [log(U(end,1)*total_cell_num)];
                y2 = [log(U(end,2)*total_cell_num)];
            else
                y1 = [U(end,1)*total_cell_num];
                y2 = [U(end,2)*total_cell_num];
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
                [u_new,v_new] = radiotherapy_no_srvn(U,LQ_param);
                [T,U] = ode45(@(t,U) stem_ODE_feedback_no_srvn(t, U, r1, r2, d, p, h, z, l, n),[sim_resume_days(i), treat_days(i+1)],[u_new v_new]);
                x = [x T(1,1) T(end,1)];
                if logQ
                    y1 = [y1 log(U(1,1)*total_cell_num) log(U(end,1)*total_cell_num)];
                    y2 = [y2 log(U(1,2)*total_cell_num) log(U(end,2)*total_cell_num)];
                else
                    y1 = [y1 U(1,1)*total_cell_num U(end,1)*total_cell_num];
                    y2 = [y2 U(1,2)*total_cell_num U(end,2)*total_cell_num];
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

                for i=1:len_samples
                    sample_times(i) = find(T < treat_start + pre_sample_times(i),1,'last');
                end
                if jj > 1
                    dose_reprog_time(jj,k,:) = U(sample_times,1);
                else
                    dose_reprog_time(jj,k,:) = U(sample_times+dose1offset,1);
                end

            end
%             if jj > 1
%                 dose_reprog_time(jj,k,:) = U(sample_times,1);
%             else
%                 dose_reprog_time(jj,k,:) = U(sample_times+dose1offset,1);
%             end
%            x = [];
%             y1 = [];
%             y2 = [];
            y3 = [];
            y4 = [];
            y5 = [];
            y6 = [];
            y7 = [];



        end
        if jj > 1
            times(jj,:) = T(sample_times);
        else
            times(jj,:) = T([sample_times(1)+dose1offset sample_times(2:end)]);
        end
        %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
    end
end
strC = arrayfun(@(x) num2str(x),C,'UniformOutput',false);
if weak_feedback_Q
    fdbk_type = 'Weak';
else
    fdbk_type = 'Strong';
end
title_prefix = [fdbk_type,' Fdbk; Relative CSC at time ~ '];
s_ylabel = '#CSC/#CSC(conventional)';
s_xlabel = 'Dose (Gy)';
s_legend_title = 'reprogramming chance';
save_prefix = ['C:\Users\jhvo9\Documents\vojh MCSB Repo\Projects\Yu\Varying Fractionation\Better Visuals\',lower(fdbk_type),'_relCSC_time_'];
xlim_vec = [min(Doses) max(Doses)];
ylim_vec = [0 3];
leg_pos = 'north';
line_types = cell(1,len_C);
for k=2:len_C
    line_types{k} = '--';
end
line_types{1} = '-';
suffix = '_alt'; % i.e. '_testN'
fig_start = 50;
for t=1:len_samples
    figX = figure(50+(t-1));
    hold on;
    for k=1:len_C
        plot(Doses,dose_reprog_time(:,k,t)/dose_reprog_time(2,k,t),'o',...
            'color', color{k},'LineStyle',line_types{k},'LineWidth',line_lengths(2));
    end
    lgdX = legend(strC,'Location',leg_pos);
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