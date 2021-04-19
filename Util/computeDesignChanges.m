function del = computeDesignChanges(ds)
[~, ld] = size(ds);
del = zeros(ld-1, 1);
dk = ds(:, 1);
for k = 1:ld-1
    dkp1 = ds(:, k+1);
    del(k) = norm(dkp1 - dk);
    dk = dkp1;
end

end