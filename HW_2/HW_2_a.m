% MRAC Simulation using MIT Rule
clc;
clear;

% Parameters
b = 2;            % Plant parameter
gamma = 1;       % Adaptation gain
K0 = 1;           % Initial controller gain
Tsim = 10;        % Simulation time
dt = 0.01;       % Time step
time = 0:dt:Tsim-dt; % Time vector

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
uc = sin(2*pi*0.2*time);         % Step reference input
dKdt = zeros(size(time));      % Adaptive rate of change of K
syms 's';
k=3;
% Simulation loop
for k = 2:length(time)
    % Compute reference model output
    tspan = [time(k-1), time(k)];
    [t, y_m] = lsim(ref_model, uc(k-1:k), tspan);
    ym(k) = y_m(end);
    
    % Compute plant output
    [t, y_p] = lsim(plant, K * (uc(k-1:k) - y(k-1)), tspan);
    y(k) = y_p(end);
    
    % Compute error
    e(k) = ym(k) - y(k);
    
    % Compute sensitivity (∂e/∂K)
    denom = (s^2 + 5*s + b*K);  % Denominator in sensitivity
    sensitivity = (b * s * (s + 5)) / denom^2;
    
    % Adaptive law (MIT Rule)
    dKdt(k) = gamma * e(k) * ilaplace(sensitivity);
    K = K + dKdt(k) * dt; % Update K
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
plot(time, dKdt, 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('K (Adaptive Gain)');
title('Adaptive Gain K');
grid on;
