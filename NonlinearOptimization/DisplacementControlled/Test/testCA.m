%% Running tests
dzGuesses = [1 0 0];
maxitss = [4 4 0];
ntests = numel(dzGuesses);
for k = 1:ntests
    dzGuess = dzGuesses(k);
    maxits = maxitss(k);
    SymmetricDC;
end

%% Test solution
nits = 201;
f1 = figure(1);
ax1 = nexttile;
f2 = figure(2);
ax2 = nexttile;
ntests = 3;
ind = 20:201;
indsol = 1:3;
load('solutions.mat');
c0 = -30.9;
colors = ['g', 'g', 'k'];
hold on;
load('solutions.mat');
for i = 1:ntests
    j = indsol(i);
    data = solutions{j};
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
    
    plot(ax1, nlqi, nfi, [sign, color], 'Displayname', name);
    hold(ax1, 'on');
    
    semilogy(ax2, ind, abs(g0i(ind)), [line, color], ...
        'Displayname', name);
    hold(ax2, 'on');
end
xlabel(ax1, 'Number of lineqs')
ylabel(ax1, 'Number of factorizations')
xlabel(ax2, 'Number of iterations')
ylabel(ax2, 'Compliance')
figure(f1), legend('Location', 'best'), figure(f2), legend('Location', 'best')