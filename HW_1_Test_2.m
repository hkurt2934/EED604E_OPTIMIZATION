% Recursive Least Square Algorithm (RLS) with Forgetting Factor
N = 1000;
forgetting_factor = 0.98;

% Desired values:
a_desired = 1; 
b_desired = 0.5;
c_desired = -0.5;

U = load('u.mat');
u = U(1).u;
Y = load("y.mat");
y = Y(1).y;

% White Gaussian noise
e = sqrt(0.5) * randn(N, 1);

% Simulate the system: y(t) = -a*y(t-1) + b*u(t-1) + e(t) + c*e(t-1)
for t = 2:N
    y(t) = -a_desired * y(t-1) + b_desired * u(t-1) + c_desired * e(t-1) + e(t);
end

theta_hat = zeros(3, 1); % Initial parameter estimates [a, b, c]
P = 100 * eye(3); % Initial covariance matrix
phi = zeros(3, 1); % Regression vector

% Recursive estimation loop with forgetting factor
theta_history_forgetting = zeros(3, N); % To store parameter estimates over time
for t = 2:N
    % Form the regression vector
    phi = [-y(t-1); u(t-1); e(t-1)];
    
    % Compute prediction error
    epsilon = y(t) - phi' * theta_hat;
    
    % Compute Kalman gain with forgetting factor
    K = P * phi / (forgetting_factor + phi' * P * phi);
    
    % Update parameter estimates
    theta_hat = theta_hat + K * epsilon;
    
    % Update covariance matrix with forgetting factor
    P = (1 / forgetting_factor) * (P - K * phi' * P);
    
    % Store parameter estimates
    theta_history_forgetting(:, t) = theta_hat;
end

% Plot the results for RLS with forgetting factor
figure;
plot(2:N, theta_history_forgetting(1, 2:N), 'r', 'LineWidth', 1.5); hold on;
plot(2:N, theta_history_forgetting(2, 2:N), 'g', 'LineWidth', 1.5);
plot(2:N, theta_history_forgetting(3, 2:N), 'b', 'LineWidth', 1.5);
yline(a_desired, '--r', 'True a');
yline(b_desired, '--g', 'True b');
yline(c_desired, '--b', 'True c');
xlabel('Time step');
ylabel('Parameter estimate');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Best');
title('Parameter Estimation using RLS with Forgetting Factor');
grid on;

