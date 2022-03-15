%% Test CG
L = 1;
N = 100;
dx = L/(N + 1);
xx = linspace(0, L, N + 2)';
xxi = xx(2:end-1);
A = -toeplitz([-2 1 zeros(1, N - 2)])/dx^2;
n = 3;
b = -(n*pi)^2*sin(n*pi*xxi);

xs = A\b;

figure();
ax = nexttile();
hold(ax, 'ON');
plot(ax, xxi, xs, 'g');

R = chol(A) + eye(size(A))*1;
R = eye(size(A));

x0 = xs + 0.1*sin(pi*xxi);
[xcg, nits] = pCG(@(x) A*x, b, R);

plot(ax, xxi, xcg, 'r');