%% Load and compute data
load solutions
refData = solutions{1};
designs = refData.stats.designs;

zold = designs(:, 1);
zold = zold/norm(zold);

coth = ones(201, 1);
for k = 2:201
    znew = designs(:, k);
    znew = znew/norm(znew);
    coth(k) = znew'*zold;
    zold = znew;
end
cothmax = 1 - 1e-3;
sinthmax = sqrt(1 - cothmax^2);
%% Plotting
f = figure;
ax = nexttile;
plot(ax, sqrt(1-coth(2:end).^2));
hold on;
% plot([1 201], [sinthmax sinthmax], 'LineStyle', '--', 'Color', [.2 .2 .2])
axis tight;
xlabel('Iteration');
ylabel('sin\theta^i', 'Interpreter', 'tex');
title('Sine of angle between previous and current design');
set(gcf, 'Position', [100 100 800 450])