function V = eNormOrth(R, K)
% Orthonormalizes the basis R with the basis K such that vi'Kvj = d_ij
[~, s] = size(R);
r1 = R(:, 1);
v1 = r1/(r1'*K*r1)^(0.5);
V = [v1];
for i = 2:s
    ri = R(:, i);
    vi = ri;
    for j = 1:i-1
        vj = V(:, j);
        vi = vi - (ri'*K*vj)*vj;
    end
    vi = vi/(vi'*K*vi)^(0.5);
    V(:, i) = vi;
end
end