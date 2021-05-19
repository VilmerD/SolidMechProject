%% Running tests
dzGuesses = [1 0 0];
maxitss = [6 6 0];
ntests = numel(dzGuesses);
for k = 1:ntests
    dzGuess = dzGuesses(k);
    maxits = maxitss(k);
    DC;
end

%% Test solution
nits = 201;
f1 = figure(1);
f2 = figure(2);
ntests = 3;

load('solutions.mat');
c0 = -30.9;
colors = ['g', 'g', 'k'];
hold on;
load('solutions.mat');
c0 = -30.9;
for i = 1:ntests
    data = solutions{end-ntests+i};
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