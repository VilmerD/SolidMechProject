f = figure;
stats = data.stats;
sampledSensitivities = stats.sampledSensitivities;
[samples, iterations] = size(sampledSensitivities);
numSens = sampledSensitivities(2:samples/2, :);
analSens = sampledSensitivities((samples/2 + 2):end, :);
relDiff = (numSens - analSens)./numSens;
ax = nexttile;
plot(ax, (1 - relDiff)');
xlabel('Iteration number')
ylabel('relative differance')
diffmax = max(abs(relDiff(:)));
e = ceil(log10(diffmax));
markminus = sprintf('1-1e%1i', e);
markplus = sprintf('1+1e%1i', e);
yticks([1 - 10^(e), 1, 1 + 10^(e)]);
yticklabels({markminus, 1, markplus})
axis([0, iterations, 1 - 10^(e), 1 + 10^(e)])