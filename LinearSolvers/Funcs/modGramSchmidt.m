function R = modGramSchmidt(R)
% Modified gram schmidt orthogonlaization
[~, s] = size(R);
for k = 1:s
    vk = R(:, k);
    rkk = norm(vk);
    qk = vk/rkk;
    for j = k + 1:s
        vj = R(:, j);
        rkj = qk'*vj;
        R(:, j) = vj - rkj*qk;
    end
end
end