%%
load solutions
refData = solutions{2};
refSensitivities = refData.stats.sensitivities;
testData = solutions{1};
testSensitivities = testData.stats.sensitivities;
nbasis = testData.nbasis;
maxits = testData.maxits;
reldiff = abs((refSensitivities - testSensitivities)./refSensitivities);

load(refData.geomfile);
%% Animate and make movie
% Animation
clear P A
f = figure();
ax = nexttile;
sForm = 'it: %03i max: %1.2e';
frames = cell2struct(cell(201, 2), {'cdata', 'colormap'}, 2);
for i = 1:201
    % Compute relative difference and plot using fill function
    if exist('P', 'var')
        setPatchData(P, 'CData', reldiff(:, i));
    else
        P = fill(ax, ex', ey', reldiff(:, i), 'EdgeColor', 'none');
        % Adjusting labels and such
        axis image;
        ax.XAxis.TickLength = [0 0];
        ax.YAxis.TickLength = [0 0];
        xticklabels({});
        yticklabels({});
        set(f, 'Position', ...
            [700; 900;
            600; 250]);
        cbar = colorbar(ax);
        caxis([0 1])
    end
    drawnow;
    
    s = sprintf(sForm, i, max(reldiff(:, i)));
    if exist('A', 'var')
        A.String = s;
    else
        A = annotation('textbox', [0.13 0.87 0.1 0.1], 'string', s, ...
            'FitBoxToText', 'on');
    end
    
    % Save frame and clear figure for next plot
    figure(f);
    frame = getframe(f);
    frames(i) = frame;
end

% Make movie
nameform = 'Sensitivities_ff%02i_nb%02i.avi';
name = sprintf(nameform, maxits, nbasis);
v = VideoWriter(name);
v.FrameRate = 6;
open(v);
writeVideo(v, frames);
close(v);