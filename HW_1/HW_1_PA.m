% Projection Algorithm
N = 1000;
alpha_pa = 0.1;
gamma_pa = 0.01;

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
phi = zeros(3, 1);

theta_store = zeros(3, N);

for t = 2:N
    phi = [-y(t-1); u(t-1); e(t-1)];
    error = y(t) - phi' * theta_hat;
    K = (gamma_pa*phi)/(alpha_pa + phi'*phi);
    theta_hat = theta_hat + K*error;
    theta_store(:, t) = theta_hat;
end

% Plot the results
figure;
plot(theta_store(1, 2:N), 'r', 'LineWidth', 2); 
hold on;
plot(theta_store(2, 2:N), 'g', 'LineWidth', 2);
plot(theta_store(3, 2:N), 'b', 'LineWidth', 2);
yline(a_desired, '--r', 'Desired a');
yline(b_desired, '--g', 'Desired b');
yline(c_desired, '--b', 'Desired c');
xlabel('Sample');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Best');
title('Projection Algorithm');
grid on;
f = gcf;
exportgraphics(f,'PA.png');