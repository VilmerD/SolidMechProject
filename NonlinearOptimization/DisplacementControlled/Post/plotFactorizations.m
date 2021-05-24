%% Plots the number of factorizations vs the number of lineqs
nits = 201;
f1 = figure;
ax1 = nexttile;
f2 = figure;
ax2 = nexttile;
imin = 180; imax = 201;
ind = imin:imax;
indsol = 1:2:17;
load('solutions.mat');
nameform = ['$f_{\\textnormal{fact.}}$: %02i, ', ...
    '$\\textnormal{n}_{\\textnormal{basis}}$: %02i'];
signs = 'sph*';  ns = length(signs);
colors = 'rbgm'; nc = length(colors);
factmax = 1;
lineqmax = 1;
for i = 1:length(indsol)
    j = indsol(i);
    data = solutions{j};
    stats = data.stats;
    g0i = stats.g0;
    ii = 200;
    fi = stats.factorizations;
    nfi = sum(fi(1:ii));
    lqi = stats.lineqs;
    nlqi = sum(lqi(1:ii));
    
    if data.maxits == 0
        color = 'k';
        line = '^--';
        sign = '^';
    else
        color = colors(ceil(i/2));
        sign = signs(mod(i-1, 2) + 1);
        line = [sign, '-'];
    end
    
    name = sprintf(nameform, data.maxits, data.nbasis);
    plot(ax1, nlqi, nfi, [sign, color], 'Displayname', name);
    hold(ax1, 'on');
    
    semilogy(ax2, ind, abs(g0i(ind)), [line, color], ...
        'Displayname', name, 'MarkerIndices', 1:10:(imax - imin));
    hold(ax2, 'on');
    factmax = max(factmax, nfi);
    lineqmax = max(lineqmax, nlqi);
end
xlabel(ax1, 'Number of lineqs')
ylabel(ax1, 'Number of factorizations')
xlabel(ax2, 'Iteration')
ylabel(ax2, 'Compliance')
% hold(ax1, 'off')
% axis(ax1, [0 lineqmax 0 factmax]);
figure(f1), leg1 = legend('Location', 'best', 'Interpreter', 'Latex');
figure(f2), leg2 = legend('Location', 'best', 'Interpreter', 'Latex');