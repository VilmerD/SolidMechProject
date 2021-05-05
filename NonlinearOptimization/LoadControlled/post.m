%% Numerical vs Analytical Derivatives
figure(1)
dc = stats.del_dc;
[nnod, ~] = size(dc);
nnodhalf = round(nnod/2, 0);
numerical_sens = dc(1:nnodhalf, :);
analytical_sens = dc(nnodhalf+1:end, :);

delta = numerical_sens - analytical_sens;
relative_delta = delta./analytical_sens; 
semilogy(abs(relative_delta'))

xlabel("Iteration number")
ylabel("Difference in derivatives")
title("Difference in analytical and numerical derivative")
hold off;

%% Number of factorizations
figure(2);
n = countNumberSolves(stats.designUpdate);
f = countFactorizations(stats.factorizations, stats.designUpdate);
coth = computeDesignChanges(stats.designs);
ymax = max(max(n), max(f));
ymin = 0;
xmax = max(length(n), length(f));
xmin = 1;
subplot(211)
axis([xmin, xmax, ymin, ymax]);
hold on;
plot(f, 'g', 'Displayname', '# of factorizations');
plot(n, 'r', 'Displayname', '# of lineqs. solved');
xlabel("Iteration number")
legend();
hold off;

subplot(212);
semilogy(coth, 'b', 'Displayname', 'Similarity in design');
axis([1, xmax, min(coth), 1]);
xlabel("Iteration number")
legend();