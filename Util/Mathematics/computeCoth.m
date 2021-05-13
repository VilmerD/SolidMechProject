function coth = computeCoth(designs)
[~, n] = size(designs);
coth = zeros(n-1, 1);
dnm1 = designs(:, 1);
dnm1 = dnm1/norm(dnm1);
for k = 1:n-1
    dnp1 = designs(:, k+1);
    dnp1 = dnp1/norm(dnp1);
    coth(k) = dnm1'*dnp1;
    dnm1 = dnp1;
end
end