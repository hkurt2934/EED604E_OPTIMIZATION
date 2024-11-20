% Extended Least Squares (ELS) Algorithm

N = 1000;
sigma2 = 0.5;

a_desired = 1; 
b_desired = 0.5;
c_desired = -0.5;

u = ones(N, 1);
y = zeros(N, 1);
E = load("e.mat");
e = E(1).e;

for t = 2:N
    y(t) = -a_desired * y(t-1) + b_desired * u(t-1) + c_desired * e(t-1) + e(t);
end
theta_hat = zeros(3, 1); % [a, b, c]
P = 100 * eye(3);
phi = zeros(3, 1);

theta_store = zeros(3, N);
for t = 2:N
    phi = [-y(t-1); u(t-1); e(t-1)];
    error = y(t) - phi' * theta_hat;
    M = (P*phi)/(1 + phi'*P*phi);
    P = P - (P*phi*phi'*P) / (1 + phi'*P*phi);
    theta_hat = theta_hat + M * error;
    theta_store(:, t) = theta_hat;
end

figure;
plot(1:iter, theta_store(1, 1:iter), 'r', 'LineWidth', 2); 
hold on;
plot(1:iter, theta_store(2, 1:iter), 'g', 'LineWidth', 2);
plot(1:iter, theta_store(3, 1:iter), 'b', 'LineWidth', 2);
yline(a_desired, '--r', 'Desired a');
yline(b_desired, '--g', 'Desired b');
yline(c_desired, '--b', 'Desired c');
xlabel('Number of Iteration');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Best');
title('ELS Algorithm');
grid on;
f = gcf;
exportgraphics(f,'ELS_1.png');
