function myeldisp2(ex, ey, ed, mark)
edx = ed(:, 1:2:end);
edy = ed(:, 2:2:end);

ex = ex + edx;
ey = ey + edy;

explot = [ex ex(:, 1)]';
eyplot = [ey ey(:, 1)]';

if nargin < 4
    mark = 'k';
end
if ischar(mark)
    plot(explot, eyplot, mark);
else
    map = [1 1 1
           0 0 0];
    fill(ex, ey, mark);
    colormap(map);
end
end