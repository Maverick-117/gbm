function [u_new,v_new, s_new] = radiotherapy(U,LQ_para, surv_vec)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u = U(end,1);
v = U(end,2);
s = U(end,3);
[a1, b1, a2, b2, c, D] = LQ_para{:};
[cont_p_a, cont_p_b, compt_mult, srvn_csc, srvn_dcc] = surv_vec{:};
% compt_mult tunes fold-change of stem cell feedback ratio based on diff.
% cell feedback ratio
u_new = u*exp(-a1*fdbk(cont_p_a, s)*D-b1*fdbk(cont_p_b*compt_mult, s)*D^2) + min(1,c*D)*v*exp(-a2*fdbk(cont_p_a, s)*D-b2*fdbk(cont_p_b, s)*D^2); % apply RT
% u*exp(-a1*D-b1*D^2) + min(1,c*D)*v
            % u*exp(-a1*D-b1*D^2) + min(1,c*D)*v*exp(-a2*D-b2*D^2);
            % min(1,c*D) is one is better due to consistency
            % established between u_new and v_new
v_new = (v - min(1,c*D)*v)*exp(-a2*fdbk(cont_p_a, s)*D-b2*fdbk(cont_p_b, s)*D^2); % apply RT; assumes death occurs on reprogrammed cells. is this biologically valid?
% max(v*exp(-a2*D-b2*D^2) - c*v*D,0)
s_new = s + srvn_csc * (u-u*exp(-a1*fdbk(cont_p_a, s)*D-b1*fdbk(cont_p_b*compt_mult, s)*D^2)) + srvn_dcc * (v-v_new);
end

function val = fdbk(control, surv)
val = 1/(1+control*surv);
end