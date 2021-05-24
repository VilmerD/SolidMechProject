%% Get the model data so the plot can be made
load(data.geomfile, 'ex', 'ey');

[nz, ndesigns] = size(data.stats.designs);
g0 = data.stats.g0;
g1 = data.stats.g1;

% Make a frame for each of the designs
f = figure(1);
subplot(122)
for k = 1:ndesigns
    ax = nexttile;
    dk = data.stats.designs(:, k);
    fill(ex', ey', dk);
    drawnow;
    axis image;
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    xticklabels({});
    yticklabels({});
    set(f, 'Position', ...
        [700; 900;
        600; 250]);
    s = sprintf('it: %03i, g0: %4.2f, g1: %1.3e', k, g0(k), g1(k));
    annotation('textbox', [0.13 0.87 0.10 0.1], 'string', s)
    frame = getframe(f);
    frames(k) = frame;
    clf(f);
    fprintf('Frame: %i\n', k)
end

% Make movie
name = 'PoorQuality2.avi';
v = VideoWriter(name);
v.FrameRate = 6;
open(v);
writeVideo(v, frames);
close(v);