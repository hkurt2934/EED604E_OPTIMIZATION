%%
a = [-0.01, 0, 0.01];
K = 1;
L_1 = 1+a(1);
L_2 = 1+a(2);
L_3 = 1+a(3);
M_1 = a(1);
M_2 = a(2);
M_3 = a(3);
num = 1;
den_1 = [K L_1 M_1];
den_2 = [K L_2 M_2];
den_3 = [K L_3 M_3];
sys_1 = tf(num,den_1);
sys_2 = tf(num,den_2);
sys_3 = tf(num,den_3);
figure(1);
step(sys_1);
xlim([0 150]);
ylim([0 400]);
hold on;
step(sys_2);
xlim([0 300]);
ylim([0 400]);
hold on;
step(sys_3);
xlim([0 300]);
ylim([0 400]);
figure(2);
sys_close_1 = feedback(sys_1,1);
sys_close_2 = feedback(sys_2,1);
sys_close_3 = feedback(sys_3,1);
step(sys_close_1);
xlim([0 10]);
ylim([0 1.2]);
hold on;
step(sys_close_2);
xlim([0 10]);
ylim([0 1.2]);
hold on;
step(sys_close_3);
xlim([0 10]);
ylim([0 1.2]);



%%
%Hope this could help you!
% Define the transfer function
num = [10 120 660 2280 4730 4600 1600];
den = [1 14 113 628 2379 6350 11115 9800 3200];
sys = tf(num, den);
% Plot the unit step response
step(sys);
% Add title and axis labels
title('Unit Step Response');
xlabel('Time (seconds)');
ylabel('Output');
% Customize the axis
xlim([0 10]); % set x-axis limits
ylim([0 1.2]); % set x-axis limits
grid on; % add grid lines
syms s;
