clc;
clear;

% Parameters
b = 2;            % Plant parameter
gamma = 10;       % Adaptation gain
K0 = 1;           % Initial controller gain
Tsim = 10;        % Simulation time
dt = 0.001;       % Time step
time = 0:dt:Tsim; % Time vector

% Reference model parameters (Gm(s) = 9 / (s^2 + 5s + 9))
Am = [1 5 9];
Bm = [9];
ref_model = tf(Bm, Am);

% Plant parameters (G(s) = b / (s(s+5)))
Ap = [1 5 0];
Bp = [b];
plant = tf(Bp, Ap);

% Initialize variables
K = K0;                        % Initial adaptive gain
e = zeros(size(time));         % Tracking error
ym = zeros(size(time));        % Reference model output
y = zeros(size(time));         % Plant output
uc = ones(size(time));         % Step reference input
dK = zeros(size(time));        % Adaptive rate of change of K

% Simulation loop
for k = 2:length(time)
    % Compute reference model output
    tspan = [time(k-1), time(k)];
    [t, y_m] = lsim(ref_model, uc(k-1:k), tspan, ym(k-1));
    ym(k) = y_m(end);
    
    % Compute plant output
    [t, y_p] = lsim(plant, K * (uc(k-1:k) - y(k-1)), tspan, y(k-1));
    y(k) = y_p(end);
    
    % Compute error
    e(k) = ym(k) - y(k);
    
    % Compute sign of error (absolute error criterion)
    sign_e = sign(e(k));
    
    % Approximate sensitivity ∂e/∂K
    % Sensitivity approximation: delta_y / delta_K
    delta_K = 0.01; % Small perturbation for numerical differentiation
    perturbed_K = K + delta_K;
    
    % Perturbed plant output
    [t, y_p_perturbed] = lsim(plant, perturbed_K * (uc(k-1:k) - y(k-1)), tspan, y(k-1));
    y_perturbed = y_p_perturbed(end);
    
    % Sensitivity approximation (numerical derivative)
    sensitivity = -(y_perturbed - y(k)) / delta_K;
    
    % Adaptive law (MIT Rule for Absolute Error)
    dK(k) = -gamma * sign_e * sensitivity;
    K = K + dK(k) * dt; % Update adaptive gain
end

% Plot results
figure;
subplot(3,1,1);
plot(time, ym, 'b', 'LineWidth', 1.5); hold on;
plot(time, y, 'r', 'LineWidth', 1.5);
legend('y_m (Reference Model)', 'y (Plant)');
xlabel('Time (s)');
ylabel('Output');
title('System Output Tracking');
grid on;

subplot(3,1,2);
plot(time, e, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Error (e)');
title('Tracking Error');
grid on;

subplot(3,1,3);
plot(time, dK, 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('K (Adaptive Gain)');
title('Adaptive Gain K');
grid on;
