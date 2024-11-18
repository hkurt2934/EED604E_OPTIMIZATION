% Parameters
N = 100; % Number of samples
sigma2 = 0.5; % Variance of noise e(t)

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

% Extended Least Squares (ELS) initialization
theta_hat = zeros(3, 1); % Initial parameter estimates [a, b, c]
phi = zeros(3, 1); % Regression vector
max_iter = 50; % Maximum number of iterations
tol = 1e-6; % Convergence tolerance

% Storage for tracking parameter estimates
theta_history = zeros(3, max_iter);

% Extended Least Squares Iterative Algorithm
for iter = 1:max_iter
    Phi = []; % Design matrix
    Y = [];   % Output vector
    
    % Construct the regression model based on current parameter estimates
    for t = 2:N
        % Estimate e(t-1) using current parameters
        if t > 2
            e_prev_est = y(t-1) - (-theta_hat(1)*y(t-2) + theta_hat(2)*u(t-2) + theta_hat(3)*e(t-2));
        else
            e_prev_est = 0; % No estimate for e(t-1) at the start
        end
        
        % Form the regression vector
        phi = [-y(t-1); u(t-1); e_prev_est];
        
        % Update the design matrix and output vector
        Phi = [Phi; phi'];
        Y = [Y; y(t)];
    end
    
    % Solve for the new parameter estimates using least squares
    theta_new = (Phi' * Phi) \ (Phi' * Y);
    
    % Check for convergence
    if norm(theta_new - theta_hat) < tol
        break;
    end
    
    % Update the parameter estimates
    theta_hat = theta_new;
    theta_history(:, iter) = theta_hat;
end

% Trim unused iterations
theta_history = theta_history(:, 1:iter);

% Plot the results
figure;
plot(1:iter, theta_history(1, :), 'r', 'LineWidth', 1.5); hold on;
plot(1:iter, theta_history(2, :), 'g', 'LineWidth', 1.5);
plot(1:iter, theta_history(3, :), 'b', 'LineWidth', 1.5);
yline(true_a, '--r', 'True a');
yline(true_b, '--g', 'True b');
yline(true_c, '--b', 'True c');
xlabel('Iteration');
ylabel('Parameter estimate');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Best');
title('Parameter Estimation using Extended Least Squares (ELS)');
grid on;
