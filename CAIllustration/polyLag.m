function l = generateLagrangian(xp, yp)
    function y = L(xx)
        lx = numel(xp);
        xxv = ones(lx - 1, 1)*xx;
        
        y = 0;
        for k = 1:lx
            xk = xp(k);
            xkv = setdiff(xp, xk);
            lk = prod((xxv - xkv'), 1)/(prod(xk - xkv));
            y = y + lk*yp(k);
        end
    end
l = @L;
end