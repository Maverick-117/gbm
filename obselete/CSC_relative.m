% other plot (interesting)
% with reprogramming 
y1 = [9.3078, 8.9355, 8.727, 8.5634, 8.5172, 8.5269, 8.6372, 8.6996, 8.771];
y2 = y1./8.727;
y3 = [5.8144 5.9473 5.95 6.001 6.016 6.0315 6.107 6.106 6.1222];
y4 = y3./5.95;
x = [1, 1.5, 2, 3, 4, 5, 6, 7, 8];

figure(1)
hold on
plot(x,y2,'color', 'r','LineWidth', 2)
plot(x,y4,'color', 'k','LineWidth', 2)
legend('reprogramming','no reprogramming')
title('\fontsize{14} Relative CSC at t = 500')
ylabel('CSC#/CSC#(conventional)')
xlabel('Dose (Gy)')
hold off



