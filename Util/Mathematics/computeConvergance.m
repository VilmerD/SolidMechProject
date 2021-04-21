function cs = computeConvergance(residuals)
lr = length(residuals);
cs = zeros(lr, 1);
cs(1) = 1;

rkm1 = residuals(1);
for k = 2:lr-1
    rk = residuals(k + 1);
    if rkm1 > 0
        cs(k) = rk/rkm1;
    end
    rk = rkm1;
end

end