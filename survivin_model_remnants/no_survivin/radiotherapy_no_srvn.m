function [u_new,v_new] = radiotherapy(U,LQ_para)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u = U(end,1);
v = U(end,2);
[a1, b1, a2, b2, c, D] = LQ_para{:};
u_new = u*exp(-a1*D-b1*D^2) + min(1,c*D)*v*exp(-a2*D-b2*D^2); % apply RT
% u*exp(-a1*D-b1*D^2) + min(1,c*D)*v
            % u*exp(-a1*D-b1*D^2) + min(1,c*D)*v*exp(-a2*D-b2*D^2);
            % min(1,c*D) is one is better due to consistency
            % established between u_new and v_new
v_new = (v - min(1,c*D)*v)*exp(-a2*D-b2*D^2); % apply RT; assumes death occurs on reprogrammed cells. is this biologically valid?
% max(v*exp(-a2*D-b2*D^2) - c*v*D,0)
end