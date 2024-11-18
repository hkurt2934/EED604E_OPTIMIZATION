clear; clc;
n = 2;
mo = 1;

X = [-1 1 -1 1];
Y = [0.44 -0.648 0.481 -0.615];

parameter_a = -1*ones(1,1000);
parameter_b = ones(1,1000);
output = [parameter_b;parameter_a];
l = 0.1;
N = 1000;
bn = zeros(n,N);
yp = [];
xn = zeros(n,1);
P = 1e5*eye;

for i=mo+1:N
    xn(:,i) = [-Y(i-1:-1:i-mo)' X(i-1:-1:i-mo)'];
    xx = inv(l)*P*xn(:,i);
    K = inv(1+xn(:,i)'*xx)*xx;
    P = (inv(l)*P)-xx'*K;
    yp(i,1) = xn(:,i)'*bn(:,i-1);
    error(i,1) = Y(i,1) - yp(i,1);
    bn(:, i) = bn(:,i-1) + K*error(i,1);
end







