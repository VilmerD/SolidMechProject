%% Get the model data so the plot can be made
load(data.geomfile, 'ex', 'ey');

[nz, ndesigns] = size(data.stats.designs);
g0 = data.stats.g0;
g1 = data.stats.g1;

% Make a frame for each of the designs
clear P A
f = figure(1);
frames = cell2struct(cell(201, 2), {'cdata', 'colormap'}, 2);
ax = nexttile;
for k = 1:ndesigns
    dk = data.stats.designs(:, k);
    if exist('P', 'var')
        setPatchData(P, 'CData', dk);
    else 
        P = fill(ex', ey', dk, 'EdgeColor', 'none');
        axis image;
        ax.XAxis.TickLength = [0 0];
        ax.YAxis.TickLength = [0 0];
        xticklabels({});
        yticklabels({});
        set(f, 'Position', ...
            [700; 900;
            600; 250]);
    end
    drawnow;
    
    s = sprintf('it: %03i, g0: %4.2f, g1: %1.3e', k, g0(k), g1(k));
    if exist('A', 'var')
        set(A, 'String', s);
    else
        A = annotation('textbox', [0.13 0.87 0.10 0.1], 'string', s);
    end
    
    frame = getframe(f);
    frames(k) = frame;
    fprintf('Frame: %i\n', k)
end

% Make movie
name = 'FineQuality3.avi';
v = VideoWriter(name);
v.FrameRate = 6;
open(v);
writeVideo(v, frames);
close(v);