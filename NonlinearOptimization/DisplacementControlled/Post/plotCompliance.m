%% 
if ~exist('data', 'var')
    load solutions.mat;
    data = solutions{end};
    clear solutions
end
g0 = data.stats.g0;
f = figure();
ax = nexttile; 
xx = 1:201;
plot(ax, xx, -g0(xx));
xlabel('Iteration'), ylabel('Compliance')
title('Compliance over iteration');
axis tight

%% 
f = figure()
ax = nexttile;
xx = 175:201;
plot(ax, xx, -g0(xx))
set(ax, 'Position', [50 50 80 80])
axis tight