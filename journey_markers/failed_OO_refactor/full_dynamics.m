function [output_no_treatment,output_treatment,output_y] = full_dynamics(prm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Defining Variables

load('fit_result_data_GBM.mat')
par = parameters;
% prm = {g, DT, D, frac_num, treat_start, acq_days_after_RT, p, h, l, z, n, 
%           ROI_radius, rho, srv_start, c,
%           total_start_frac, F, saveQ, logQ, cont_p_a, cont_p_b,
%           compt_mult, srvn_csc, srvn_dcc}
g = prm.g; DT = prm.DT; p = prm.p; h = prm.h; l = prm.l; z = prm.z; n = prm.n;
srv_start = prm.srv_start; logQ = prm.logQ; c = prm.c;


total_cell_num = 4/3*pi()*prm.ROI_radius^3*prm.rho;
F = par(g,1);
a1 = par(g,2);  
b1 = par(g,3); 
a2 = par(g,4); 
b2 = par(g,5); 
% dedifferentiation rate
r1 = log(2)/DT; % growth rate of CSC
r2 = log(2)/DT; % growth rate of DCC
d  = log(2)/DT; % death rate of DCC
sc_start = prm.total_start_frac*F;
tc_start = total_start_frac-sc_start;
% Defining treatment days and ODE simulation start time after each fraction
weeks = floor(prm.frac_num/5);
total_days = prm.frac_num + 2*weeks;
acq_end = prm.treat_start + total_days + prm.acq_days_after_RT - 1;
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
x = [T(end,1)];
T_r = T; Ta_r = Ta;

surv_vec = {prm.cont_p_a, prm.cont_p_b, prm.compt_mult, prm.srvn_csc, prm.srvn_dcc};

if logQ
    y1 = [log10(U(end,1)*total_cell_num)];
    y2 = [log10(U(end,2)*total_cell_num)];
    ys = [log10(U(end,3)*total_cell_num)];
    U_r = log10(U(1:end-1,:));
    Ua_r = log10(Ua(1:end-1,:));
else
    y1 = [U(end,1)*total_cell_num];
    y2 = [U(end,2)*total_cell_num];
    ys = [U(end,3)*total_cell_num];
    U_r = U(1:end-1,:);
    Ua_r = Ua(1:end-1,:);
end
y3 = [U(end,1)];
y4 = [U(end,2)];
S = U(:,1)+U(:,2);
y5 = [U(end,1)./S(end)*100];
y6 = [p./(1+l.*(U(end,2)).^n)];
y7 = [r1./(1+h.*(U(end,2)).^z)];


for i = 1:length(sim_resume_days)
%%%%%%% with stem cell %%%%%%%%%
    LQ_param = {a1, b1, a2, b2, c, D};

    [u_new,v_new, s_new] = radiotherapy(U,LQ_param, surv_vec);

    [T,U] = ode45(@(t,U) stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig),[sim_resume_days(i), treat_days(i+1)],[u_new v_new s_new]);
    x = [x T(1,1) T(end,1)];
    T_r = [T_r; T]; 
    if logQ
        y1 = [y1 log10(U(1,1)*total_cell_num) log10(U(end,1)*total_cell_num)];
        y2 = [y2 log10(U(1,2)*total_cell_num) log10(U(end,2)*total_cell_num)];
        ys = [ys log10(U(1,3)*total_cell_num) log10(U(end,3)*total_cell_num)];
        U_r = [U_r; log10(U(1:end-1,:))];
    else
        y1 = [y1 U(1,1)*total_cell_num U(end,1)*total_cell_num];
        y2 = [y2 U(1,2)*total_cell_num U(end,2)*total_cell_num];
        ys = [ys U(1,3)*total_cell_num U(emd,3)*total_cell_num];
        U_r = [U_r; U(1:end-1,:)];
    end
    y3 = [y3 U(1,1) U(end,1)];
    y4 = [y4 U(1,2) U(end,2)]; 
    SS = U(:,1)+U(:,2);
    y5 = [y5 U(1,1)./SS(1)*100 U(end,1)./SS(end)*100];
    y6 = [y6 p./(1+l.*(U(1,2)).^n) p./(1+l.*(U(end,2)).^n)];
    y7 = [y7 r1./(1+h.*(U(1,2)).^z) r1./(1+h.*(U(end,2)).^z)];

    xx = x(1:end-1);
    yy = y6(1:end-1);
    yyy = y7(1:end-1);
    finalU = U;
end
output_no_treatment = {Tb, Ub};
output_treatment = {[Ta_r; T_r], [Ua_r; U_r]};
        
output_y = {x, y1, y2, y3, y4, SS, y5, y6, y7, finalU, ys};

end

