keepvars = {'relPlotQ','absTimeQ','relPlotQVals','absTimeQVals'};
delvars = setdiff(who,keepvars);
relPlotQVals = [true false];
absTimeQVals = [true false];
for relPlotQ=relPlotQVals
    for absTimeQ=absTimeQVals
        %% Basic 
        clc;
        keepvars = {'relPlotQ','absTimeQ','relPlotQVals','absTimeQVals'};
        delvars = setdiff(who,keepvars);
        clear(delvars{:},'delvars');
        close all
        load('fit_result_data_GBM.mat')

        sig_vec = [-1 0.1 1 10] * sig;
        len_sig = length(sig_vec);
        %% Defining Variables

        relPlotQVals = [true false];
        absTimeQVals = [true false];

        %%%%%              SWITCHBOARD START              %%%%%
        %relPlotQ = false; % true means normalize by conventional; false means just plot the value
        %absTimeQ = false; % true means plot everything according to an absolute time; false means plot N days after end of radiotherapy
        saveQ = true; % true means plots are saved; false means plots are not saved
        logQ = true; % [DEPRECATED; holdover from time-series plot.]
        weak_feedback_Q = true; % true means use a very low feedback gain on self-renewal probability; false means use a very high feedback gain on same constant
        srvQ = true; % true means that survivin is at play for any sim; false means that survivin is barred from any sim
        %%%%%               SWITCHBOARD END               %%%%%

        % Radiotherapy model
        par = parameters;
        DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
        Doses = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % dose fraction sizes in Gy
        Frac = [60 30 20 15 12 10 8 7 6 6 5 5 4 4 4 3 3 3 3 3]; 
        len_D = length(Doses);
        %sample_times = [1 5 10 20 30]; 
        pre_pre_sample_times = [0 70 100 300 600 900 1200 1500];
        %dose1offset = 0;
        line_lengths = [0.8 2.3];
        treat_start = 100; % treatment start time relative to start of simulation
        acq_days_after_RT = 2500; % simulation length after last fraction
        len_samples = length(pre_pre_sample_times);

        % ODE model

        p = .505; % self renewal probability 

        h = 10^(5); z = 1; % control the strenth of inhibitory signal
        %low = -3; high = -3; c_l = 1*(high-low+1);
        %C = [0 5.196*10.^(linspace(low,high,c_l))]; % ramp up by factors of 2 later
        C = [5.196*10^(-3)];
        % [0 5.196*2.^(linspace(low,high,c_l))];
        len_C = length(C);
        times = zeros(len_D,len_samples);
        offset = 600;

        mu_bar = C; chi = 10^4;


        if weak_feedback_Q
            l = 10^(-7); % weak feedback
        else
            l = 10^(3); % strong feedback 
        end
        n = 1; 



        % dedifferentiation rate
        r1 = log(2)/DT; % growth rate of CSC
        r2 = log(2)/DT; % growth rate of DCC
        d  = log(2)/DT; % death rate of DCC

        % Initial Condition
        total_start_frac = 0.2/64;
        ROI_radius = 1; % Radius of the simulation region of intrest
        rho = 10^9; % density of cells in the region of interest (cells/cm^3)
        total_cell_num = 4/3*pi()*ROI_radius^3*rho;
        total_start_cell_num = total_cell_num*total_start_frac;

        zeta_mult1 = 0.1;
        zeta_mult2 = 1;
        if srvQ
            %srvn_csc = 3.6; srvn_dcc = 0.05;
            srvn_zeta = [3.6 * zeta_mult1, 0.05 * zeta_mult2]; %; 3.6, 0.5; 3.6, 5];
            % model is more sensitive to this than to cont_p
        else
            %srvn_csc = 0; srvn_dcc = 0;
            srvn_zeta = [0, 0];
        end


        cont_p_a = 10^4;
        cont_p = [1 10] * cont_p_a;
        %cont_p = [0 0; 1 0.1; 1 1; 1 10; 1 100] * cont_p_a;
        % the simulation seems to be sort of sensitive to control parameters
        %cont_p = [1 10] * cont_p_a;
        %cont_p_a = 10^4; % \gamma_{\alpha_V}
        %cont_p_b = [0.1 10]; % \gamma_{\beta_V}
        compt_mult = 100; % \zeta_U; in reality, this should be higher than 1 since stem cells 
                        % tend to more difficult to kill, so a higher denominator is more attractive
        len_surv_cont = size(cont_p); 
        len_surv_cont = len_surv_cont(1);
        len_zeta = size(srvn_zeta); len_zeta = len_zeta(1);
        dose_reprog_time_cont = zeros(3,len_D,len_C,len_samples,len_zeta,len_sig,len_surv_cont);
        N_c = len_sig;
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

                        for ss = 1:len_sig
                            sig = sig_vec(ss);
                            if sig <= 0
                                surv_vec = {cont_p(ii,1), cont_p(ii,2), compt_mult, 0, 0}; %assuming these control parameters are constantd
                                sig = 0;
                            else
                                surv_vec = {cont_p(ii,1), cont_p(ii,2), compt_mult, srvn_zeta(gg,1), srvn_zeta(gg,2)}; %assuming these control parameters are constantd
                            end
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
                                [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[0, treat_days(1)],[sc_start tc_start srv_start], options);
                                [Ta,Ua] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[0, treat_days(1)],[sc_start tc_start srv_start], options);
                                % Without treatment 
                                [Tb,Ub] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[0, acq_end],[sc_start tc_start srv_start], options);
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
        %                             if sig > 10^6
        %                                 s_new = 0;
        %                             end
                                    [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig, mu_bar, chi),[sim_resume_days(i), treat_days(i+1)],[u_new v_new s_new]);
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

                                if absTimeQ
                                    % measurement at absolute time
                                    pre_sample_times = treat_start + Frac(1) + pre_pre_sample_times; % 
                                    pre_sample_times(1) = pre_sample_times(1) + 30;
                                else
                                    % measurement at time after last radiotherapy
                                    % session
                                    pre_sample_times = sim_resume_days(end) + pre_pre_sample_times; 
                                end
                                sample_times = 1:len_samples;
                                for i=1:len_samples
                                    sample_times(i) = find(T <=  pre_sample_times(i),1,'last'); % 0 vs ;
                                end
                                %sample_times = find(T <= sim_resume_days(end) + offset, 1, 'last');
                                dose_reprog_time_cont(1,jj,k,:,gg,ss,ii) = U(sample_times,1);%U(sample_times,2) * mu_bar * chi * U(sample_times,3)./(1+chi * U(sample_times,3));%
                                dose_reprog_time_cont(2,jj,k,:,gg,ss,ii) = U(sample_times,2);
                                dose_reprog_time_cont(3,jj,k,:,gg,ss,ii) = U(sample_times,3);
        %                         if jj > 1
        %                             dose_reprog_time_cont(jj,k,:,gg,ss,ii) = mu_bar * chi * U(sample_times,3)/(1+chi * U(sample_times,3));
        %                         else
        %                             dose_reprog_time_cont(jj,k,:,gg,ss,ii) = U(sample_times+dose1offset,3);
        %                         end
                    %            x = [];
                    %             y1 = [];
                    %             y2 = [];
                                y3 = [];
                                y4 = [];
                                y5 = [];
                                y6 = [];
                                y7 = [];
                                times(jj,:) = T(sample_times);
        %                         if jj > 1
        %                             times(jj,:) = T(sample_times);
        %                         else
        %                             times(jj,:) = T([sample_times(1)+dose1offset sample_times(2:end)]); %same as before but with offset at first time to let times line up across dosage.
        %                         end


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
        s_ylabel = 'csc';%'dediff * V';%/dediff (conventional)';

        s_xlabel = 'Dose (Gy)';
        s_legend_title = 'Survivin Decay Rate:';
        save_prefix = ['C:\Users\jhvo9\Google Drive (vojh1@uci.edu)\a PhD Projects\GBM Modeling\Survivin Models\Model 2 - Dedifferentiation Term\optimality_studies\' num2str(srvn_zeta(1)) '_' num2str(srvn_zeta(2)) '\'];%'C:\Users\jhvo9\Documents\vojh MCSB Repo\Projects\Yu\Varying Fractionation\Better Visuals\',lower(fdbk_type), '_relCSC_time_'];
        xlim_vec = [min(Doses) max(Doses)];
        ylim_vec = [-Inf Inf];%[-Inf Inf];
        leg_pos = 'north';
        line_types = cell(1,len_surv_cont);
        for k=2:len_surv_cont
            line_types{k} = '--';
        end
        line_types{1} = '-';
        if relPlotQ
            pre_suffix = '_rel_';
        else
            pre_suffix = '_abs_';
        end
        if absTimeQ
            post_suffix = '_abs_';
        else
            post_suffix = '_rel_';
        end
        suffix = [pre_suffix s_ylabel post_suffix 'time'];%'_abs_total_abs_time'; % i.e. '_testN'
        if relPlotQ
            s_ylabel = [s_ylabel  '/' s_ylabel '(conventional)'];
        end
        per_ratio = zeros(len_D,len_C,len_samples,len_surv_cont);
        fig_start = 102;
        for gg=1:len_zeta
            for t=1:len_samples
                title_prefix = {[fdbk_type,' Fdbk; (' num2str(cont_p(ii,1)) ',' num2str(cont_p(ii,2)) ',' num2str(compt_mult) ');'] ['(' num2str(srvn_zeta(gg,1)) ',' num2str(srvn_zeta(gg,2)) '); ' 'Total Pop at ']};% ]}; ' at end of treatment + ' num2str(offset) 'days'
                figX = figure(fig_start*gg+(t-1));
                hold on;
                for ss=1:len_sig
                    if sig_vec(ss) <= 0
                        disp_str = "no survivin";
                    else
                        disp_str = num2str(sig_vec(ss));
                    end
                    for k=1:len_C
                        for ii=1:len_surv_cont
                            % when relQ = true, we normalize by the conventional
                            % treatment-associated value.
                            % when relQ = false, we simply do nothing more to it.
                            ordinate = zeros(length(Doses),3);
                            for d=1:len_D
                                ordinate(d,:) = (dose_reprog_time_cont(1,d,k,t,gg,ss,ii));
                                
                                % worth measuring:
                                % csc: (dose_reprog_time_cont(1,d,k,t,gg,ss,ii)
                                % SF: exp(-a1./(1 + cont_p(ii,1)*compt_mult* dose_reprog_time_cont(3,d,k,t,gg,ss,ii))*Doses(d)-b1./(1 + cont_p(ii,2)*compt_mult* dose_reprog_time_cont(3,d,k,t,gg,ss,ii))*Doses(d)^2);
                                % CSC renewal: (2*p./(1+l*(dose_reprog_time_cont(2,d,k,t,gg,ss,ii))^n)-1)*r1./(1+h*(dose_reprog_time_cont(2,d,k,t,gg,ss,ii))^z)*(dose_reprog_time_cont(1,d,k,t,gg,ss,ii));
                                % CSC derivative: stem_ODE_feedback(0, dose_reprog_time_cont(:,d,k,t,gg,ss,ii), r1, r2, d, p, h, z, l, n, sig, mu_bar, chi)';
                                % total: (dose_reprog_time_cont(1,d,k,t,gg,ss,ii) + dose_reprog_time_cont(2,d,k,t,gg,ss,ii));
                                % csc-to-dcc: dose_reprog_time_cont(1,d,k,t,gg,ss,ii)./dose_reprog_time_cont(2,d,k,t,gg,ss,ii);

                                % not worth measuring:

                                % r1./(1+h.*(dose_reprog_time_cont(2,d,k,t,gg,ss,ii)).^n);
                                %mu_bar * chi * dose_reprog_time_cont(3,d,k,t,gg,ss,ii)./(1 + chi * dose_reprog_time_cont(3,d,k,t,gg,ss,ii));
                            end

                            if relPlotQ
                                ordinate = ordinate./ordinate(2,:);
                            end
                            ptemp = plot(Doses, ordinate(:,2),'o',... % dose_reprog_time_cont(:,k,t,gg,ss,ii)  * total_cell_num
                                    'color', color{ss},'LineStyle',line_types{ii},'LineWidth',line_lengths(2),...
                                    'DisplayName',disp_str,'MarkerSize',3*ss);
                            % 
                            per_ratio(:,k,t,gg,ii) = max(0,1-dose_reprog_time_cont(1,:,k,t,gg,ss,ii)/dose_reprog_time_cont(1,2,k,t,gg,ss,ii));
                            ptemp.Color(4) = 0.6;
        %                 pMarkers = ptemp.MarkerHandle; 
        %                 pMarkers.MarkerEdgeColor = uint8(255*[color{ss}'; 0.5]);
                        end
                    end
                end
                lgdX = legend('Location',leg_pos);
                title(lgdX,s_legend_title)
                %[title_prefix,num2str(round(mean(times(:,t)),ceil(log10(mean(times(:,t)))),'significant')),' days']
                if absTimeQ
                    title_prefix{2} = [title_prefix{2} 'time ~ ' num2str(round(mean(times(:,t)),ceil(log10(mean(times(:,t)))),'significant')),' days'];
                else
                    title_prefix{2} = [title_prefix{2} num2str(pre_pre_sample_times(t)) ' days post-radiotherapy end'];
                end
                title(title_prefix)
                ylabel(s_ylabel)
                xlabel(s_xlabel)
                xlim(xlim_vec)
                ylim(ylim_vec)
                hold off;                                            
                if saveQ
                    save_label = [save_prefix num2str(sample_times(t)) '_' num2str(srvn_zeta(gg,1)) '_' num2str(srvn_zeta(gg,2))];
                    savefig([ save_label suffix '.fig'])
                    saveas(gcf,[ save_label suffix '.png'])
                end
            end
        end
    end
end
        per_ratio(per_ratio == 0) = nan;
                   
