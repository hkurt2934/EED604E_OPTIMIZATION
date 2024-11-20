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
phi = zeros(3, 1);
max_iter = 50;
tolerance = 1e-6;
theta_store = zeros(3, N);
for iter = 1:max_iter
    
    Phi = zeros(N, 3);
    Y = zeros(N, 1);
    for t = 3:N
        e_prev_est = y(t-1) - (-theta_hat(1)*y(t-2) + theta_hat(2)*u(t-2) + theta_hat(3)*e(t-2));
        phi = [-y(t-1), u(t-1), e_prev_est];
        Phi(t-1, :) = phi;
        Y(t-1) = y(t);
    end
    
    theta_new = (Phi' * Phi) \ (Phi' * Y);
    if norm(theta_new - theta_hat) < tolerance
        break;
    end
    
    theta_hat = theta_new;
    theta_store(:, iter) = theta_hat;
end

theta_store = theta_store(:, 1:iter);

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
exportgraphics(f,'ELS.png');
