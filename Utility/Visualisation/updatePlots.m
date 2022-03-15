function updatePlots(g0, g1, alph, j, lines)
% UPDATEPLOTS upadtes the plots with the new information
%

% Plotting goal function values etc
g0_line = lines(1);
g0_line.YData = [g0_line.YData g0];
g0_line.XData = [g0_line.XData j];

% Plotting constraint
for k = 1:numel(g1)
    g1k_line = lines(k + 1);
    g1k_line.YData = [g1k_line.YData g1(k)];
    g1k_line.XData = [g1k_line.XData j];
end

% Plotting design
desn_line = lines(end);
desn_line.YData = [desn_line.YData alph];
desn_line.XData = [desn_line.XData j];
drawnow;
end