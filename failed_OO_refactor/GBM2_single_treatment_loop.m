%% Basic 
clc;
clear;
load('fit_result_data_GBM.mat')

%% Defining Variables
% Radiotherapy model
par = parameters;%[.052 .01 .071 .197 .203];%   
DT =  3.9; % doubling time in days in order labeled by varible "cell_lines"
Doses = [6]; % dose fraction sizes in Gy
Frac = [10]; 
treat_start = 100; % treatment start time relative to start of simulation
acq_days_after_RT = 400; % simulation length after last fraction

% ODE model
color = {'k','r','b'};
p = .505; % self renewal probability 
x =[]; y1 = []; y2 = [];
h = 10^(5); z = 1; % control the strength of inhibitory signal
weak_feedback_Q = false;
if weak_feedback_Q
    l = 10^(-7); % weak feedback
    fdbk_type = 'Weak';
else
    l = 10^(3); % strong feedback 
    fdbk_type = 'Strong';
end
n = 1; 
pwr = 3;
C = [0 5.196*10^(-pwr)];
% dedifferentiation rate
r1 = log(2)/DT; % growth rate of CSC
r2 = log(2)/DT; % growth rate of DCC
d  = log(2)/DT; % death rate of DCC

% Initial Condition
total_start_frac = 0.2/64; % Victoria: 0.2/64; Nayeon: 0.0005/64
ROI_radius = 1; % Radius of the simulation region of intrest
rho = 10^9; % density of cells in the region of interest (cells/cm^3)
total_cell_num = 4/3*pi()*ROI_radius^3*rho;
total_start_cell_num = total_cell_num*total_start_frac;
suffix = '_alt';
saveQ = false; logQ = true; srvQ = true;
srv_start = 0; %initializing survivin amount
cont_p_a = 1; cont_p_b = 0.1; compt_mult = 2; 

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
        ppre_prm = {'g' g; 'DT' DT; 'D' D; 'frac_num' frac_num;
          'treat_start' treat_start; 'acq_days_after_RT' acq_days_after_RT;
          'p' p; 'h' h; 'l' l; 'z' z; 'n' n; 'total_start_frac' total_start_frac; 
          'ROI_radius' ROI_radius; 'rho' rho; 'srv_start' srv_start; 'c' C(k);
          'saveQ' saveQ; 'logQ' logQ; 'cont_p_a' cont_p_a; 'cont_p_b' cont_p_b;
           'compt_mult' compt_mult; 'srvn_csc' srvn_csc; 'srvn_dcc' srvn_dcc};
        pre_prm = ppre_prm'; prm = struct(pre_prm{:});
        % where full_dynamics.m applies        
        [output_no_treatment,output_treatment,output_y] = full_dynamics(prm);
        Tb = output_no_treatment{1}; Ub = output_no_treatment{2};
        T = output_treatment{1}; U = output_treatment{2};
        x = output_y{1}; y1 = output_y{2}; y2 = output_y{3}; y3 = output_y{4};
        y4 = output_y{5}; SS = output_y{6}; y5 = output_y{7}; y6 = output_y{8};
        y7 = output_y{9}; U_final = output_y{10}; ys = output_y{11};
        
        figure(1)
        hold on
        if logQ
            plot(Tb(:,1),log10(Ub(:,1)*total_cell_num),'g','LineStyle','--','LineWidth', 2) % stem cell
            plot(Tb(:,1),log10(Ub(:,2)*total_cell_num),'g','LineWidth', 2) 
            ylabel('\fontsize{12} Survival cell number (log)')
        else
            plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', 2) % stem cell
            plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', 2) 
            ylabel('\fontsize{12} Survival cell number')
        end
        if weak_feedback_Q
            title({'\fontsize{12} with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{12} with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
        
        xlabel('\fontsize{12} Time (Days)')
        axis([0 acq_days_after_RT -inf inf])
        hold off 
        
        figure(100)
        hold on
        if logQ
            plot(Tb(:,1),log10((Ub(:,1)+Ub(:,2))*total_cell_num),'g','LineStyle','-.','LineWidth', 2)            
            plot(x(1:end-1),log10(10.^(y1(1:end-1))+10.^(y2(1:end-1))),'color', color{k},'LineStyle','-.','LineWidth', 2)
            plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel('\fontsize{12} Survival cell number (log)')
        else
            plot(T(:,1),(Ub(:,1)+Ub(:,2))*total_cell_num,'g','LineStyle','-.','LineWidth', 2)
            ylabel('\fontsize{12} Survival cell number')
        end
        if weak_feedback_Q
            title({'\fontsize{12} Total size with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{12} Total size with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
        
        xlabel('\fontsize{12} Time (Days)')
        axis([0 acq_days_after_RT -inf inf])

        if logQ
            plot(x(1:end-1),log10(10.^(y1(1:end-1))+10.^(y2(1:end-1))),'color', color{k},'LineStyle','-.','LineWidth', 2)
            plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
        else
            plot(x(1:end-1),y1(1:end-1)+y2(1:end-1),'color', color{k},'LineStyle','-.','LineWidth', 2)
            plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
        end
            
        if logQ
            plot(T(:,1),log10((U(:,1)+U(:,2))*total_cell_num),'color', color{k},'LineStyle','-.','LineWidth', 2) % stem cell
            plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
         else
            plot(T(:,1),(U(:,1)+U(:,2))*total_cell_num,'color', color{k},'LineStyle','-.','LineWidth', 2) % stem cell
            plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        hold off 
        
        figure(10)
        hold on
        if logQ
            plot(Tb(:,1),log10(Ub(:,1)*total_cell_num),'g','LineStyle','-.','LineWidth', 2) % stem cell
            plot(Tb(:,1),log10(Ub(:,2)*total_cell_num),'g','LineWidth', 2) 
            plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            plot(Ta(:,1),log10(Ua(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel({'Survival cell'; 'number (log)'})
        else
            plot(Tb(:,1),Ub(:,1)*total_cell_num,'g','LineStyle','--','LineWidth', 2) % stem cell
            plot(Tb(:,1),Ub(:,2)*total_cell_num,'g','LineWidth', 2) 
            plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            plot(Ta(:,1),Ua(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
            ylabel({'Survival cell'; 'number'})
        end
        xlabel('Time (Days)')
        axis([0 acq_days_after_RT -inf inf])
        if logQ
            plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)           
        else
            plot(T(:,1),U(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),U(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', 2)
        plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', 2)
        hold off
     
        figure(2)
        hold on
        plot(Tb(:,1),Ub(:,1)./(Ub(:,1)+Ub(:,2))*100,'g','LineWidth', 2) 
        plot(x(1:end-1),y5(1:end-1),'color', color{k},'LineWidth', 2)
        plot(T(:,1),U(:,1)./(U(:,1)+U(:,2))*100,'color', color{k},'LineWidth', 2)
        if weak_feedback_Q
            title({'\fontsize{12} with feedback on m_u,m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{12} with feedback on m_u,m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
         
%        
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} CSC percentage (%)')
        axis([0 acq_days_after_RT .1 200])
        hold off
        
        figure(11)
        hold on
        plot(Tb(:,1), Ub(:,1)./(Ub(:,1)+Ub(:,2))*100,'g','LineWidth', 2) 
        plot(x(1:end-1),y5(1:end-1),'color', color{k},'LineWidth', 2)
        plot(T(:,1),U(:,1)./(U(:,1)+U(:,2))*100,'color', color{k},'LineWidth', 2) 
        axis([0 acq_days_after_RT .1 200])
        xlabel('Time (Days)')
        ylabel({'CSC'; 'percentage (%)'})
        hold off
        
        figure(4)
        hold on
        axis([0 acq_days_after_RT -inf inf])
        if weak_feedback_Q
            title({'\fontsize{14} Self Renewal Probability';'\fontsize{12} feedback on m_u, m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{14} Self Renewal Probability';'\fontsize{12} feedback on m_u, m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
         
%        
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} p value')
        plot(Tb(:,1),p./(1+l.*(Ub(:,2)).^n),'g','LineWidth', 2)
        plot(Ta(:,1),p./(1+l.*(Ua(:,2)).^n),'color', color{k},'LineWidth', 2) 
        plot(x(1:end-1),y6(1:end-1),'color',color{k},'LineWidth', 2)
        
        tt = T(:,1);
        pp = p./(1+l.*(U(:,2)).^n);
        plot([xx(end),tt(1)],[yy(end),pp(1)],'color',color{k},'LineWidth', 2)
        plot(T(:,1),p./(1+l.*(U(:,2)).^n),'color',color{k},'LineWidth', 2)
        hold off
        
        figure(5)
        hold on
        axis([0 acq_days_after_RT 0.025 inf])
        if weak_feedback_Q
            title({'\fontsize{14} Division Rate';'\fontsize{12} feedback on m_u, m_v and p (weak)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        else
            title({'\fontsize{14} Division Rate';'\fontsize{12} feedback on m_u, m_v and p (strong)'; ['Dose:',num2str(Doses(1)),' Gy, Days:',num2str(Frac(1))]})
        end
        xlabel('\fontsize{12} Time (Days)')
        ylabel('\fontsize{12} m_u = m_v value')
        plot(Tb(:,1),r1./(1+h.*(Ub(:,2)).^z),'g','LineWidth', 2)
        plot(x(1:end-1),y7(1:end-1),'color',color{k},'LineWidth', 2)
        rr = r1./(1+h.*(U(:,2)).^z);
        plot([xx(end),T(1)],[yyy(end),rr(1)],'color',color{k},'LineWidth', 2)
        plot(T(:,1),r1./(1+h.*(U(:,2)).^z),'color',color{k},'LineWidth', 2)
        hold off
         
            
        figure(1)
        hold on
        if logQ
            plot(Ta(:,1),log10(Ua(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),log10(Ua(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            %plot(Ta(:,1),log10(Ua(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
            plot(T,log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
        else
            plot(Ta(:,1),Ua(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(Ta(:,1),Ua(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            plot(T,U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
            %plot(Ta(:,1),Ua(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        plot(x(1:end-1),y1(1:end-1),'color', color{k},'LineStyle','--','LineWidth', 2)
        plot(x(1:end-1),y2(1:end-1),'color', color{k},'LineWidth', 2)
        %plot(x(1:end-1),ys(1:end-1),'color', color{k},'LineStyle',':','LineWidth', 2)
        
        if logQ
            plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
            plot(T(:,1),log10(U(:,1)*total_cell_num),'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),log10(U(:,2)*total_cell_num),'color', color{k},'LineWidth', 2) 
            plot(T(:,1),log10(U(:,3)*total_cell_num),'color', color{k},'LineStyle',':','LineWidth', 2)
        else
            plot(T(:,1),U(:,1)*total_cell_num,'color', color{k},'LineStyle','--','LineWidth', 2) % stem cell
            plot(T(:,1),U(:,2)*total_cell_num,'color', color{k},'LineWidth', 2) 
            plot(T(:,1),U(:,3)*total_cell_num,'color', color{k},'LineStyle',':','LineWidth', 2)
        end
        hold off 

        




        
    end
end
    %clearvars -except sc_start tc_start p color tUV_all tUVb_all treat_all sim_resume_all TCP_all TCPb_all par DT Doses Frac treat_start ROI_radius total_cell_num total_start_frac total_start_cell_num rho acq_days_after_RT cell_lines data_new par_ab parameters F a b a1 a2 b1 b2 g k TCPb_frac TCP_frac BED BEDb TCP_thre folder
        
fig1 = figure(1);
fig2 = figure(10);
[h_m h_i]=inset(fig1,fig2);
set(h_i,'xtick',100:20:160,'xlim',[99,160])

fig3 = figure(2);
fig4 = figure(11);
[h_m h_i]=inset(fig3,fig4);
set(h_i,'xtick',100:20:160,'xlim',[99,160],'yscale','log')
set(h_m,'yscale','log')

cab 3 4 5 6 100
fig1 = figure(3);
fig2 = figure(4);
fig3 = figure(5);
fig4 = figure(6);

if weak_feedback_Q
    fdbk_type = 'Weak';
else
    fdbk_type = 'Strong';
end
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