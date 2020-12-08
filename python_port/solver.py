# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 14:44:06 2020

@author: jhvo9
"""

class solver:
    def __init__(self,params):
        self.params = params;
    
    def radiotherapy_schedule(self,idx):
        g = idx[0];
        frac_num = self.Frac(g);
        sc_start = self.total_start_frac*self.F;
        # Defining treatment days and ODE simulation start time after each fraction
        weeks = int(frac_num/5);
        total_days = frac_num + 2*weeks;
        acq_end = self.treat_start + total_days + self.acq_days_after_RT - 1;
        A = [1, 1, 1, 1, 1, 0, 0] * (weeks + 1)#repmat([1 1 1 1 1 0 0], 1, weeks+1);
        A_new = A[range(total_days)];
        idx_flt = list(map(lambda x: x == 1, [0, 1, 0, 1]));
        treat_days = list(map(lambda x: x+self.treat_start-1,A_new[idx_flt]));
        sim_resume_days = treat_days+10/(60*24); # ODE simulation resumes 10 minutes after RT
        treat_days.append(acq_end);
    
"""
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
"""