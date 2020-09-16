function dU = stem_ODE_feedback(t, U, r1, r2, d, p, h, z, l, n, sig)
dU = zeros(3, 1);
% U(1) := stem cell; U(2) := diff. cell; U(3) := survivin
% V = U(2)/(U(1)+U(2));

% Yu eq (feedback on r1, r2 and p)
dU(1) = (2*p/(1+l*U(2)^n)-1)*r1/(1+h*U(2)^z)*U(1); % stem cell
dU(2) = 2*(1-p/(1+l*U(2)^n))*r1/(1+h*U(2)^z)*U(1)+r2/(1+h*U(2)^z)*U(2)-d/(1+h*U(2)^z)*U(2); %differentiated cell
dU(3) = -sig*U(3);

% Yu eq (feedback on r1, r2 and p)
% dU(1) = (2*p/(1+l*U(2)^n)-1)*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0); % stem cell
% dU(2) = 2*(1-p/(1+l*U(2)^n))*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0)+r2/(1+h*U(2)^z)*U(2)*cutoff(U(2),0)-d/(1+h*U(2)^z)*U(2)*cutoff(U(2),0); %differentiated cell


% Kim eq (feedback on r1, r2 and p)
% dU(1) = (2*p/(1+l*U(2)^n)-1)*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0); % stem cell
% dU(2) = 2*(1-p/(1+l*U(2)^n))*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0)+r2/(1+h*U(2)^z)*U(2)*cutoff(U(2),0)-d/5 *U(2)*cutoff(U(2),0); %differentiated cell


% % feedback r1, r2 and p (as fraction)   
% dU(1) = (2*p/(1+l*V^n)-1)*r1/(1+h*V^z)*U(1)*cutoff(U(1),0); 
% dU(2) = 2*(1-p/(1+l*V^n))*r1/(1+h*V^z)*U(1)*cutoff(U(1),0)+r2/(1+h*V^z)*U(2)*cutoff(U(2),0)-d*U(2)*cutoff(U(2),0);

% % feedback r1, r2 and p (as fraction but only for p)
% dU(1) = (2*p/(1+l*V^n)-1)*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0);
% dU(2) = 2*(1-p/(1+l*V^n))*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0)+r2/(1+h*U(2)^z)*U(2)*cutoff(U(2),0)-d*U(2)*cutoff(U(2),0);


% no feedback 
% % dU(1) = (2*p-1)*r1*U(1)*cutoff(U(1),0); % stem cell
% dU(2) = 2*(1-p)*r1*U(1)*cutoff(U(1),0)+r2*U(2)*cutoff(U(2),0)-d*U(2)*cutoff(U(2),0);
% feedback on r1, r2
% dU(1) = (2*p-1)*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0); % stem cell
% dU(2) = 2*(1-p)*r1/(1+h*U(2)^z)*U(1)*cutoff(U(1),0)+r2/(1+h*U(2)^z)*U(2)*cutoff(U(2),0)-d*U(2)*cutoff(U(2),0);
% feedback on p
% dU(1) = (2*p/(1+hh*U(2)^zz)-1)*r1*U(1)*cutoff(U(1),0); % stem cell
% dU(2) = 2*(1-p/(1+hh*U(2)^zz))*r1*U(1)*cutoff(U(1),0)+r2*U(2)*cutoff(U(2),0)-d*U(2)*cutoff(U(2),0); 
% k(w) exist
% dU(1) = (2*p-1)*r1/(1+h*U(2)^z)*max(1-(U(1)+U(2))^4,0)*U(1)*cutoff(U(1),0); % stem cell
% dU(2) = 2*(1-p)*r1/(1+h*U(2)^z)*max(1-(U(1)+U(2))^4,0)*U(1)*cutoff(U(1),0)+r2/(1+h*U(2)^z)*max(1-(U(1)+U(2))^4,0)*U(2)*cutoff(U(2),0)-d*U(2)*cutoff(U(2),0); 

function G = cutoff(value, cut_num)
if value > cut_num
    G = 1;
else
    G = 0;
end

