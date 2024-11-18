% Parameters
N = 100; % Number of samples
sigma2 = 0.5; % Variance of noise e(t)

% True parameters (for simulation)
true_a = 1; 
true_b = 0.5;
true_c = -0.5;

% Generate noise e(t)
e = sqrt(sigma2) * randn(N, 1); % White Gaussian noise with N(0, sigma2)

% Initialize system output y(t) and input u(t)
y = zeros(N, 1);
u = ones(N, 1);

% Simulate the system with feedback: u(t) = -0.2 * y(t)
for t = 2:N
    u(t) = -0.2 * y(t); % Feedback input
    y(t) = -true_a * y(t-1) + true_b * u(t-1) + true_c * e(t-1) + e(t);
end

% Recursive Least Squares (RLS) initialization
theta_hat = zeros(3, 1); % Initial parameter estimates [a, b, c]
P = 1e6 * eye(3); % Initial covariance matrix
phi = zeros(3, 1); % Regression vector

% Recursive estimation loop
theta_history = zeros(3, N); % To store parameter estimates over time
for t = 2:N
    % Form the regression vector
    phi = [-y(t-1); u(t-1); e(t-1)];
    
    % Compute prediction error
    epsilon = y(t) - phi' * theta_hat;
    
    % Compute Kalman gain
    K = P * phi / (phi' * P * phi + sigma2);
    
    % Update parameter estimates
    theta_hat = theta_hat + K * epsilon;
    
    % Update covariance matrix
    P = P - K * phi' * P;
    
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
title('Parameter Estimation using RLS with Feedback');
grid on;
