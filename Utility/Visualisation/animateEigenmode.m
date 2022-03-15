function frames = animateEigenmode(model, x, P, name)
% Filter
xf = model.density_filter.forward(x);
xtol = 0.1;
ifilled = xf > xtol;
xffilled = xf(ifilled);

% Amplitude
mode_amp = 10e-3;
Pa = P*mode_amp/max(abs(P));
Pae = Pa(model.edof(ifilled, 2:end));
Pax = Pae(:, 1:2:end);
Pay = Pae(:, 2:2:end);

% Find a good scaling of the image
ex = model.ex(ifilled, :);
ey = model.ey(ifilled, :);
w = max(ex(:)) - min(ex(:));
h = max(ey(:)) - min(ey(:)) + 2*mode_amp;
q = w/h;

% Resolution in height
rh = 200;

% Make figure and axis
fig = figure('Name', 'Untitled');
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
axis(ax, [0 w -mode_amp h+mode_amp]);
xticks(ax, {});
yticks(ax, {});
hold(ax, 'ON');

% Create one frame per design
T = 5;
frame_rate = 20;
N = T*frame_rate;
tt = linspace(0, 1, N);
a = sin(2*pi*tt);

frames = cell2struct(cell(N, 2), {'cdata', 'colormap'}, 2);
for k = 1:N
    ak = a(k);
    if k == 1
        Pa = fill(ax, ex', ey', xffilled, 'EdgeColor', 'none');
    else
        exk = ex + Pax*ak;
        eyk = ey + Pay*ak;
        setPatchData(Pa, {'XData', 'YData'}, {exk', eyk'});
    end
    drawnow;
    
    frames(k) = getframe(fig);
    fprintf('%03i/%03i (%02i%%)\n', k, N, round(k/N*100, 0));
end
repframes = repmat(frames, 10, 1);

% Make movie
v = VideoWriter(name);
v.FrameRate = frame_rate;
open(v);
writeVideo(v, repframes);
close(v);
end

function setPatchData(P, fields, data)
%SETPATCHDATA sets the CData property of patches in P to data
%   Speeds the animation as new patch objects don't need to be created each
%   time the image is rendered

for i = 1:numel(fields)
    fi = fields{i};
    di = data{i};
    for j = 1:numel(P)
        set(P(j), fi, di(:, j));
    end
end
end