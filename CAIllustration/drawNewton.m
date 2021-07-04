function drawNewton(ax, x, y, k)
hold(ax, 'ON');
for i = 1:numel(k)
    xi = x(i);
    xip1 = x(i + 1);
    xxi = linspace(xi, xip1);
    yi = y(i);
    li = yi + (xxi - xi) * k(i);
    plot(ax, xxi, li, 'k--');
end
end