%%
load solutions
data1 = solutions{end-1};
sens1 = data1.stats.sensitivities;
designs = data1.stats.designs;

data2 = solutions{end};
sens2 = data2.stats.sensitivities;

load(data1.geomfile);
%% Animate and make movie
% Animation
clear P A
f = figure();
ax = nexttile;
sForm = 'it: %03i absmax: %1.3f';
frames = cell2struct(cell(201, 2), {'cdata', 'colormap'}, 2);
for i = 1:201
    % Compute relative difference and plot using fill function
    reldiff = (sens1(:, i) - sens2(:, i))./sens1(:, i);
    if exist('P', 'var')
        setPatchData(P, 'EdgeColor', [abs(reldiff), zeros(nelm, 2)]);
        setPatchData(P, 'CData', designs(:, i))
    else
        P = fill(ax, ex', ey', designs(:, i));
        % Adjusting labels and such
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
    
    s = sprintf(sForm, i, max(abs(reldiff)));
    if exist('A', 'var')
        A.String = s;
    else
        A = annotation('textbox', [0.13 0.87 0.1 0.1], 'string', s);
    end
    
    % Save frame and clear figure for next plot
    figure(f);
    frame = getframe(f);
    frames(i) = frame;
end

% Make movie
name = 'Sensitivities.avi';
v = VideoWriter(name);
v.FrameRate = 6;
open(v);
writeVideo(v, frames);
close(v);