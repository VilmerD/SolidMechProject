function fig = plotDesign(model, xk, varargin)
% plotDesign plots the design xk using the fill function, for the geometry
% in model.

narginextra = numel(varargin);
if mod(narginextra, 2)
    errorStruct.message = 'Variable arguments mismatch in number';
    error(errorStruct);
end

validargs = {'Filter', 'Name', 'Visible', 'Colormap', 'Invert', ...
    'RemoveVoid'};
defaultargs = {true, 'untitled', true, 'parula', false, false};
nvalidargs = numel(validargs);
for k = 1:nvalidargs
    arg = validargs{k};
    locarg = strcmpi(varargin, arg);
    if sum(locarg)
        optargs.(arg) = varargin{find(locarg) + 1};
    else
        optargs.(arg) = defaultargs{k};
    end
end

fig = figure('Name', optargs.Name, 'Visible', optargs.Visible); 
set(fig, 'Menubar', 'None', 'Toolbar', 'None');

% Scale the image properly
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

% Selecting args
if optargs.Filter
    xf = model.density_filter.forward(xk);
else
    xf = xk;
end

xtol = 0.1;
ex = model.ex;
ey = model.ey;
if optargs.RemoveVoid
    ifilled = find(xf >= xtol);
    xf = xf(ifilled);
    ex = ex(ifilled, :);
    ey = ey(ifilled, :);
end

if optargs.Invert
    xf = 1 - xf;
end

fill(ax, ex', ey', xf, 'Linestyle', 'none');
colormap(optargs.Colormap)
axis(ax, 'tight');
xticks(ax, {});
yticks(ax, {});
end