% Parameters
N = 100; % Number of samples
sigma2 = 0.5; % Variance of noise e(t)
mu = 0.01; % Step size for LMS (choose appropriately for convergence)

% True parameters (for simulation)
true_a = 1; 
true_b = 0.5;
true_c = -0.5;

% Generate input u(t) and noise e(t)
u = randn(N, 1); % Random input signal
e = sqrt(sigma2) * randn(N, 1); % White Gaussian noise with N(0, sigma2)

% Initialize system output y(t)
y = zeros(N, 1);

% Simulate the system: y(t) + a*y(t-1) = b*u(t-1) + e(t) + c*e(t-1)
for t = 2:N
    y(t) = -true_a * y(t-1) + true_b * u(t-1) + true_c * e(t-1) + e(t);
end

% LMS initialization
theta_hat = zeros(3, 1); % Initial parameter estimates [a, b, c]
theta_history = zeros(3, N); % To store parameter estimates over time

% LMS algorithm loop
for t = 2:N
    % Form the regression vector
    phi = [-y(t-1); u(t-1); e(t-1)];
    
    % Compute prediction
    y_pred = phi' * theta_hat;
    
    % Compute prediction error
    epsilon = y(t) - y_pred;
    
    % Update parameter estimates using LMS rule
    theta_hat = theta_hat + mu * epsilon * phi;
    
    % Store parameter estimates
    theta_history(:, t) = theta_hat;
end

% Plot the results
figure;
plot(2:N, theta_history(1, 2:N), 'r', 'LineWidth', 1.5); hold on;
plot(2:N, theta_history(2, 2:N), 'g', 'LineWidth', 1.5);
plot(2:N, theta_history(3, 2:N), 'b', 'LineWidth', 1.5);
yline(true_a, '--r', 'True a');
yline(true_b, '--g', 'True b');
yline(true_c, '--b', 'True c');
xlabel('Time step');
ylabel('Parameter estimate');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Best');
title('Parameter Estimation using Least Mean Squares (LMS)');
grid on;
