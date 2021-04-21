res = stats.residuals;
% cs = computeConvergance(stats.residuals);
ifact = find(stats.factorizations > 0);
lif = length(ifact);
Du = stats.designUpdate;

figure();
semilogy(res, '-r');
hold on;
semilogy([1; 1]*ifact, [1e5; 1e-9]*ones(1, lif), '--k');
semilogy(Du, ones(length(Du), 1), 'bs');

%%
z = stats.z;
semilogy(abs(z(1:8, :))')
hold on;
semilogy([1; 1]*ifact, [1; 1e-9]*ones(1, lif), '--k');