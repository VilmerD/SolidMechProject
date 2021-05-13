%% Plot the objective function
figure(f1)
g0 = stats.g0;  g1 = stats.g1;
ax = subplot(2, 1, 1);
plot(ax, g0, 'Displayname', 'g0')
ymax = max(g0)*1.1; ymin = min(g0)*1.1;
axis(ax, [0, nits, ymin, ymax]);
xlabel('Iteration');
legend();

ax = subplot(2, 1, 2);
plot(ax, g1, 'Displayname', 'g1')
ymax = max(g1)+1e-2; ymin = min(g1)*1.3;
axis(ax, [0, nits, ymin, ymax]);
xlabel('Iteration');
legend();

sgtitle('Function value and kkt condition')
%% Number of factorizations
figure(f2);
n = stats.lineqs;
f = stats.factorizations;
coth = computeCoth(stats.designs);
ymax = max(max(n), max(f));
ymin = 0;
xmax = nits;
xmin = 1;
subplot(211)
axis([xmin, xmax, ymin, ymax]);
hold on;
plot(f, 'g', 'Displayname', '# of factorizations');
plot(n, 'r', 'Displayname', '# of lineqs. solved');
xlabel("Iteration")
legend();
hold off;

subplot(212);
semilogy(coth, 'b', 'Displayname', 'Similarity in design');
axis([1, xmax, min(coth), 1]);
xlabel("Iteration")
legend();

%% DZ plot
figure(f3)
designs = stats.designs;
[nelm, nd] = size(designs);
dz = zeros(nd-1, 1);
for k = 1:nd-1
    dz(k) = norm(designs(:, k+1) - designs(:, k));
end
semilogy(abs(dz));