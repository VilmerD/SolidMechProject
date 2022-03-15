function lines = setupPlots(g00, g10, figname)
% SETUPPLOTS prepares the plots for the optimization
% 

if nargin < 3
    figname = 'untitled_optimization';
end
fig_goal = figure('name', figname);
set(fig_goal, 'Position', [4100 400 500 600]);
tiledlayout(3, 1);

% Plot goal function
ax_goal = nexttile(1);
axis(ax_goal, 'tight');
hold(ax_goal, 'ON');

lines = plot(ax_goal, 0, g00, 'Displayname', 'g_{0}');

ylabel(ax_goal, 'g_0');
legend('Location', 'SouthWest')

% Plot constraints
ax_const = nexttile(2);
axis(ax_const, 'tight');
hold(ax_const, 'ON');

for k = 1:numel(g10)
    name = sprintf('g_{1%i}', k);
    g1k_line = plot(ax_const, 0, g10(k), 'Displayname', name);
    lines(k + 1) = g1k_line;
end

ylabel(ax_const, 'g_1');
legend('Location', 'SouthWest');

% Plot design
ax_desn = nexttile(3);
axis(ax_desn, 'tight');
set(ax_desn, 'YScale', 'log');
hold(ax_desn, 'ON');

lines(end + 1) = semilogy(ax_desn, 0, 1);

xlabel(ax_desn, 'Iteration');
ylabel(ax_desn, '\alpha', 'Fontsize', 16);
end
