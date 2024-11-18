% Parameters
N = 100; % Number of samples
sigma2 = 0.5; % Variance of noise e(t)

% True parameters (for simulation)
true_a = 1; 
true_b = 0.5;
true_c = -0.5;

% Generate input u(t) and noise e(t)
u = randn(N, 1); % Random input signal
u = ones(N, 1);
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

% Iterative solution for ELS
max_iter = 50; % Maximum number of iterations
tol = 1e-6; % Convergence tolerance

theta_history = zeros(3, max_iter); % To track parameter updates
for iter = 1:max_iter
    % Form the design matrix and output vector
    Phi = zeros(N-1, 3); % Regression matrix
    Y = zeros(N-1, 1); % Output vector
    for t = 2:N
        % Use previous estimates of noise terms
        e_prev_est = y(t-1) - (-theta_hat(1)*y(t-2) + theta_hat(2)*u(t-2) + theta_hat(3)*e(t-2));
        phi = [-y(t-1), u(t-1), e_prev_est];
        Phi(t-1, :) = phi;
        Y(t-1) = y(t);
    end
    
    % Solve the least squares problem
    theta_new = (Phi' * Phi) \ (Phi' * Y);
    
    % Convergence check
    if norm(theta_new - theta_hat) < tol
        break;
    end
    
    % Update the parameters
    theta_hat = theta_new;
    theta_history(:, iter) = theta_hat;
end

% Trim unused iterations
theta_history = theta_history(:, 1:iter);

% Plot convergence of parameters
figure;
plot(1:iter, theta_history(1, 1:iter), 'r', 'LineWidth', 1.5); hold on;
plot(1:iter, theta_history(2, 1:iter), 'g', 'LineWidth', 1.5);
plot(1:iter, theta_history(3, 1:iter), 'b', 'LineWidth', 1.5);
yline(true_a, '--r', 'True a');
yline(true_b, '--g', 'True b');
yline(true_c, '--b', 'True c');
xlabel('Iteration');
ylabel('Parameter estimate');
legend('a estimate', 'b estimate', 'c estimate', 'Location', 'Best');
title('Parameter Estimation using Extended Least Squares (ELS)');
grid on;
