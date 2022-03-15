function frames = animateIteration(model, X, name)
% Find a good scaling of the image
ex = model.ex;
ey = model.ey;
w = max(ex(:)) - min(ex(:));
h = max(ey(:)) - min(ey(:));
q = w/h;

% Resolution in height
rh = 200;
model.thresholding.b = model.thresholding.b0;

% Make figure and axis
fig = figure('Name', name);
w = max(model.ex(:));
h = max(model.ey(:));
hres = 150;
wres = ceil(w*hres/h);
set(fig, 'Position', [50 50 wres hres]);
hmargin = 0.05;
hn = 1 - 2*hmargin;
wmargin = h*hmargin/w;
wn = 1 - 2*wmargin;
ax = axes(fig, 'Position', [wmargin hmargin wn hn]);
xticklabels(ax, {});
yticklabels(ax, {});
hold(ax, 'ON');

% Create one frame per design
nsteps = size(X, 2);
frames = cell2struct(cell(nsteps, 2), {'cdata', 'colormap'}, 2);
for k = 1:nsteps
    model.thresholding.update(k);
    xk = model.density_filter.forward(X(:, k));
    if k == 1
        P = fill(ax, ex', ey', 1 - xk, 'EdgeColor', 'none');
        colormap('gray');
    else
       setPatchData(P, 'CData', 1 - xk); 
    end
    drawnow;
    
    frames(k) = getframe(fig);
    fprintf('%03i/%03i (%02i%%)\n', k, nsteps, round(k/nsteps*100, 0)); 
end

% Make movie
v = VideoWriter(['Results/Movies/', name]);
v.FrameRate = 6;
open(v);
writeVideo(v, frames);
close(v);
end

function setPatchData(P, feild, data)
%SETPATCHDATA sets the CData property of patches in P to data
%   Speeds the animation as new patch objects don't need to be created each
%   time the image is rendered
np = size(P, 1);
nd = size(data, 1);
if np ~= nd
    errorStruct.message = ...
        'Mismatch in number of patches and number of data points';
    error(errorStruct);
end
for i = 1:np
    set(P(i), feild, data(i, :));
end
end