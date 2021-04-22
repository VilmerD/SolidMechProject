function coth = computeDesignChanges(ds)
[~, ld] = size(ds);
coth = zeros(ld-1, 1);

dk = ds(:, 1);
ndk = norm(dk);

for k = 1:ld-1
    dkp1 = ds(:, k+1);
    ndkp1 = norm(dkp1);
    coth(k) = (dkp1/ndkp1)'*(dk/(ndk));
    
    dk = dkp1;
    ndk = ndkp1;
end

end