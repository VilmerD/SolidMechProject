function [K, M, Ms, T] = heatFlowMats()
% Numpoints etc
pts1d = [-sqrt(3/5), 0, sqrt(3/5)];
wts1d = [5 8 5]/9;

[Y, X] = meshgrid(pts1d);
pts = [X(:) Y(:)];

[wy, wx] = meshgrid(wts1d);
wts = reshape(wy.*wx, [], 1);

% T
N = @(x, y) 1/4*[(x - 1)*(y - 1)
                 (x + 1)*(y - 1)
                 (x + 1)*(y + 1)
                 (x - 1)*(y + 1)]';
T = zeros(4, 1);
for k = 1:size(pts, 1)
    p = pts(k, :);
    w = wts(k);
    T = T + N(p(1), p(2))'*w;
end

% M
NtN = @(x, y) N(x, y)'*N(x, y);

M = zeros(4, 4);
for k = 1:size(pts, 1)
    p = pts(k, :);
    w = wts(k);
    M = M + NtN(p(1), p(2))*w;
end

% M_surf
Ntemp = @(t) NtN(t, -1);
i = 1;
Ms = zeros(4, 1);
for k = 1:3
    p = pts1d(k);
    w = wts1d(k);
    ik = Ntemp(p)*w;
    Ms = Ms + ik(ik ~= 0);
end

% K
B = @(x, y) 1/4*[(y - 1) (y - 1) (y + 1) (y + 1)
                 (x - 1) (x + 1) (x + 1) (x - 1)];
BtB = @(x, y) B(x, y)'*B(x, y);

K = zeros(4, 4);
for k = 1:size(pts, 1)
    p = pts(k, :);
    w = wts(k);
    K = K + BtB(p(1), p(2))*w;
end
end