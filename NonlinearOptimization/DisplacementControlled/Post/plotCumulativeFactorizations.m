%% Plots the comulative number of factorizations
xx = 1:201;
index = [3 4];
step = 1;
f = figure;
signs = 'sph*';  ns = length(signs);
colors = 'rbgm'; nc = length(colors);
nameform = 'nb: %02i';
for i = 1:numel(index)
    data = solutions{index(i)};
    if data.maxits == 0
        color = 'k';
        line = '^--';
        sign = '^';
    else
        color = colors(ceil(i/step));
        sign = signs(mod(i-1, step) + 1);
        line = [sign, '-'];
    end
    name = sprintf(nameform, data.nbasis);
    plot(xx, cumsum(data.stats.factorizations), [color, line], ...
        'DisplayName', name, 'MarkerIndices', 1:10:201)
    xlabel('Iteration'), ylabel('Cumulative number of factorizations')
    hold on;
end
title('Cumulative number of factorizations for different parameter choices')
axis tight
legend('Location', 'Northwest');
set(gcf, 'Position', [100 100 800 450])