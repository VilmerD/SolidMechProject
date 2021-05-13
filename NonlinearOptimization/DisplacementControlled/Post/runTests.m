% Using dzGuess
% CA approximation
dzGuess = 1;
maxits = 4;
DC;

% No CA
maxits = 0;
DC;

% Not using dzGuess
% CA approximation
dzGuess = 0;
maxits = 4;
DC;

% No CA
maxits = 0;
DC;

%%
dzGuess_v = [1 1 0 0];
maxits_v = [4 0 4 0];
for k = 1:4
    dzGuess = dzGuess_v(k);
    maxits = maxits_v(k);
    DC;
end

%%
%% Test solution
figure(1);
hold on;
load('solutions.mat');
c0 = -12.42;
for i = 1:4
    data = solutions{end-i+1};
    stats = data.stats;
    g0i = stats.g0;
    ii = find(g0i < c0, 1, 'first');
    fi = stats.factorizations;
    nfi = sum(fi(1:ii));
    lqi = stats.lineqs;
    nlqi = sum(lqi(1:ii));
    name = sprintf('dzGuess: %i CA :%i', dzGuess_v(i), maxits_v(i) > 0);
    sign = '';
    if dzGuess_v(i) == 1
        sign = ['g', sign];
    else
        sign = ['r', sign];
    end
    if maxits_v(i) == 0
        sign = [sign, '^'];
    else
        sign = [sign, 's'];
    end
    plot(nlqi, nfi, sign, 'Displayname', name);
    hold on
end
xlabel('Number of iterations')
ylabel('Number of factorizations')
legend();