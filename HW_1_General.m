% Parameters
N = 1000;                % Number of samples
step_length = N / num_steps; % Length of each step
noise_std = 0.3;        % Standard deviation of noise

% Step function input data
u = ones(N, 1);         % Initialize input vector

% System parameters (true weights)
true_weights = [0.1; 0.5]; % Parameters for system dynamics
true_bias = 0.2;            % Bias term

% System output data (with fit to linear system dynamics)
y = zeros(N, 1);            % Initialize output vector
for k = 3:N
    % Example output model: y[k] = w1*y[k-1] + w2*u[k-2] + bias + noise
    y(k) = true_weights(1) * y(k-1) + true_weights(2) * u(k-2) + true_bias + noise_std * randn(); % Add Gaussian noise
end


%%
