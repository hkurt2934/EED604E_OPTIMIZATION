% Creating Error
N = 1000;
e = randn(N, 1);
e = abs(0.5 * e / max(abs(e)));
save('e.mat',"e");