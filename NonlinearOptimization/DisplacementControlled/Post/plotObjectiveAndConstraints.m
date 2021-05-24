%% Plot
xmin = 100; xmax = 201;
xax = xmin:xmax;
g0 = data.stats.g0;
g1 = data.stats.g1;
yyaxis left
plot(xax, abs(g0(xax)))
xlabel('Iteration number'), ylabel('Compliance')
axis tight

yyaxis right
plot(xax, g1(xax))
ylabel('Volume constraint')