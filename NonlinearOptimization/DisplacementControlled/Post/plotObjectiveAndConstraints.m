%% Plot
f = figure();
ax = nexttile;
xmin = 1; xmax = 201;
xx = xmin:xmax;
g0 = data.stats.g0;
g1 = data.stats.g1;
yyaxis left
plot(ax, xx, abs(g0(xx)))
xlabel('Iteration number'), ylabel('Compliance')
axis tight

yyaxis right
plot(ax, xx, g1(xx))
ylabel('Volume constraint')