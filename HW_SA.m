% Parameters
N = 1000; % Number of samples
sigma2 = 0.5; % Variance of noise e(t)
gamma = 0.1; % Step size for SA (learning rate)

% True parameters (for simulation)
true_a = 0.5; 
true_b = 1.2;
true_c = 0.8;

% Generate input u(t) and noise e(t)
u = randn(N, 1); % Random input signal
e = sqrt(sigma2) * randn(N, 1); % White Gaussian noise with N(0, sigma2)

% Initialize system output y(t)
y = zeros(N, 1);

% Simulate the system: y(t) + a*y(t-1) = b*u(t-1) + e(t) + c*e(t-1)
for t = 2:N
    y(t) = -true_a * y(t-1) + true_b * u(t-1) + true_c * e(t-1) + e(t);
end

% Stochastic Approximation (SA) initialization
theta_hat = zeros(3, 1); % Initial parameter estimates [a, b, c]
theta_history = zeros(3, N); % To store parameter estimates over time

% SA algorithm loop
for t = 2:N
    % Form the regression vector
    phi = [-y(t-1); u(t-1); e(t-1)];
    
    % Compute prediction error
    epsilon = y(t) - phi' * theta_hat;
    
    % Update parameter estimates using SA rule
    theta_hat = theta_hat + gamma * epsilon * phi;
    
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
title('Parameter Estimation using Stochastic Approximation (SA)');
grid on;
