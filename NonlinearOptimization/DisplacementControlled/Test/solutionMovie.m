%% Plays a movie of the solution
% Load the solution data
load solutions.mat;
solutionNumber = 11;
data = solutions{solutionNumber};
designs = data.stats.designs;

%% Get the model data so the plot can be made
load(data.geomfile, 'ex', 'ey');

[nz, ndesigns] = size(designs);
g0 = data.stats.g0;
g1 = data.stats.g1;

% Make a frame for each of the designs
f = figure(1);
subplot(122)
for k = 1:ndesigns
    ax = nexttile;
    dk = designs(:, k);
    fill(ex', ey', dk);
    drawnow;
    axis image;
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    xticklabels({});
    yticklabels({});
    set(f, 'Position', ...
        [700; 900;
        1200; 250]);
    s = sprintf('it: %03i, g0: %4.2f, g1: %1.3e', k, g0(k), g1(k));
    annotation('textbox', [0.13 0.87 0.10 0.1], 'string', s)
    frame = getframe(f);
    frames(k) = frame;
    clf(f);
    k
end

save('Movie', 'frames')
%% Play film
name = 'solutionMovie3.avi';
v = VideoWriter(name);
v.FrameRate = 6;
open(v);
writeVideo(v, frames);
close(v);