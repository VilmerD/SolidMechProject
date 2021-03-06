%% Plot
xp = [0 0.25 0.5 0.75 1];
y1 = [0, 0.35, 0.65, 0.85, 1.15];
y2 = [0, 0.65, 1.05, 1.2, 1.50];

L1 = polyLag(xp, y1);
L2 = polyLag(xp, y2);

xx = linspace(0, 1, 1000);
f = figure();
ax = nexttile;
plot(ax, xx, L1(xx), 'r', xx, L2(xx), 'r');
hold(ax, 'ON');

%% Solution
t = [1/4, 2/4, 3/4, 1];
x0 = zeros(numel(t), 1);
for i = 1:numel(t)
    ti = t(i);
    [x, k] = cheapNewton(@(x) L1(x) - ti, x0(i));
    plot(ax, [0, 1], [ti, ti], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
    plot(ax, x(end), L1(x(end)), 'k*')
    drawNewton(ax, x, L1(x), k);
    x0(i+1) = x(end);
end

ti = t(4);
[x, k] = cheapNewton(@(x) L2(x) - ti, x0(end));
plot(ax, [0, 1], [ti, ti], 'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
plot(ax, x(end), L2(x(end)), 'k*')
drawNewton(ax, x, L2(x), k);

xlabel('$u$', 'Interpreter', 'Latex');
xticks([0 1])
xticklabels([0 1])

ax.TickLabelInterpreter = 'Latex';
ylabel('$f$', 'Interpreter', 'Latex');
yticks(t)
yticklabels({'$f_1$', '$f_2$', '$f_3$', '$f_4$'})

text(0.70, 1.35, '$r(u, f, z_2) = 0$', 'Interpreter', 'Latex')
text(0.75, 0.80, '$r(u, f, z_1) = 0$', 'Interpreter', 'Latex')
text(0.92, 1.1, '$k_1^4$', 'Interpreter', 'Latex')
%%
f = figure;
ax = nexttile;
text(ax, 0.20, 0.20, '$u^{1, *}_1$', 'Interpreter', 'Latex');
text(ax, 0.40, 0.20, '$u^{2, *}_1$', 'Interpreter', 'Latex');
text(ax, 0.60, 0.20, '$u^{3, *}_1$', 'Interpreter', 'Latex');
text(ax, 0.80, 0.20, '$u^{4, *}_1$', 'Interpreter', 'Latex');
text(ax, 1.00, 0.20, '$k^{4, *}_1$', 'Interpreter', 'Latex');

text(ax, 0.20, 0.40, '$u^{4, 0}_2$', 'Interpreter', 'Latex');
text(ax, 0.40, 0.40, '$u^{4, 1}_2$', 'Interpreter', 'Latex');
text(ax, 0.60, 0.40, '$u^{4, *}_2$', 'Interpreter', 'Latex');
text(ax, 0.80, 0.40, '$k^{4, 0}_2$', 'Interpreter', 'Latex');

text(ax, 0.20, 0.60, '$CA$', 'Interpreter', 'Latex');

text(ax, 0.20, 0.80, '$f^1$', 'Interpreter', 'Latex');
text(ax, 0.40, 0.80, '$f^2$', 'Interpreter', 'Latex');
text(ax, 0.60, 0.80, '$f^3$', 'Interpreter', 'Latex');
text(ax, 0.80, 0.80, '$f^4$', 'Interpreter', 'Latex');
axis(ax, [0, 1.5, 0, 1]);