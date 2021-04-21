function N = bar3gs(ec, ed, ep)
% Computes Greens strain and the corresponding normal force
dx0 = ec(:, 1) - ec(:, 2);        % Reference vector
du = ed(1:3) - ed(4:6);         % Displacement

stretch = stretch1D(du, dx0);
sg = stress1D(ep(1), stretch);      % Internal force
N = ep(2) * sg;                     % Normal force