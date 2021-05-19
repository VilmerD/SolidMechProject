%% Running tests
dzGuess_v = [1 0 0];
maxits_v = [6 6 0];
for k = 1:4
    dzGuess = dzGuess_v(k);
    maxits = maxits_v(k);
    DC;
end

%%
%% Test solution
nits = 201;
f1 = figure(1);
f2 = figure(2);
load('solutions.mat');
c0 = -30;
colors = ['g', 'g', 'k'];
for i = 1:4
    data = solutions{end-4+i};
    stats = data.stats;
    g0i = stats.g0;
    ii = find(g0i < c0, 1, 'first');
    fi = stats.factorizations;
    nfi = sum(fi(1:ii));
    lqi = stats.lineqs;
    nlqi = sum(lqi(1:ii));
    name = sprintf('dzGuess: %i CA :%i', data.dzGuess, data.maxits);
    
    color = colors(i);
    if mod(i, 2) == 1
        line = '--';
        sign = 's';
    else
        line = '-';
        sign = '^';
    end
    
    ax1 = gca(f1);
    plot(ax1, nlqi, nfi, [sign, color], 'Displayname', name);
    hold(ax1, 'on');
    
    ax2 = gca(f2);
    semilogy(ax2, 50:201, abs(g0i(50:end)), [line, color], ...
        'Displayname', name);
    hold(ax2, 'on');
end
xlabel(ax1, 'Number of iterations')
ylabel(ax1, 'Number of factorizations')
xlabel(ax2, 'Number of iterations')
ylabel(ax2, 'Compliance')
legend();