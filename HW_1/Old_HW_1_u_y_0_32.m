% Parameters
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
    u(t) = -0.32 * y(t-1);
    y(t) = -a_desired * y(t-1) + b_desired * u(t-1) + c_desired * e(t-1) + e(t);
end

theta_hat = zeros(3, 1); % [a, b, c]
P = 100 * eye(3);
phi = zeros(3, 1);

theta_store = zeros(3, N);
for t = 2:N
    phi = [-y(t-1); u(t-1); e(t-1)];
    epsilon = y(t) - phi' * theta_hat;
    K = P * phi / (phi' * P * phi + sigma2);
    theta_hat = theta_hat + K * epsilon;
    P = P - K * phi' * P;
    theta_store(:, t) = theta_hat;
end

figure;
plot(theta_store(1, 2:N), 'r', 'LineWidth', 2);
hold on;
plot(theta_store(2, 2:N), 'g', 'LineWidth', 2);
plot(theta_store(3, 2:N), 'b', 'LineWidth', 2);
yline(a_desired, '--r', 'Desired a');
yline(b_desired, '--g', 'Desired b');
yline(c_desired, '--b', 'Desired c');
xlabel('Sample');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Southeast');
title('RLS Algorithm with 0.32*y feedback');
grid on;
