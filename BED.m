clc;
clear;
load('fit_result_data_GBM.mat');
lenD = length(Doses);
BED1 = 1:lenD;
BED2 = 1:lenD;
for i=1:lenD
    %dos = Doses(i);
    %fr = Frac(i);
    %total = fr + 0*2*floor(fr/5);
    BED1(i) = BEDcalc(Frac(i),Doses(i),a1,b1); %fr * dos * (1 + dos/(a1/b1)); 
    BED2(i) = BEDcalc(Frac(i),Doses(i),a2,b2);%fr * dos * (1 + dos/(a2/b2)); 
end
BEDavg = [mean(BED1) mean(BED2)];


% 
% classdef BED
%     properties
%         dose
%         frac_num
%         alpha
%         beta
%     end
%     methods 
%         function obj = BED(args)
%             obj.args(1) = dose;
%             obj.args(2) = frac_num;
%             obj.args(3) = alpha;
%             obj.args(4) = beta;
%         end
%         
%     end
% end