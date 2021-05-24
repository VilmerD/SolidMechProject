%% Plot g1, g0
% Load the solution data
load solutions.mat;
solutionNumber = 10;
data = solutions{solutionNumber};
designs = data.stats.designs;

%% Plot
xmin = 50; xmax = 201;
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