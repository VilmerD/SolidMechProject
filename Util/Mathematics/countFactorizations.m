function f = countFactorizations(fac, DU)
updates = length(DU) - 1;
f = zeros(updates, 1);

for k = 1:updates
    f(k) = sum(fac(DU(k):DU(k + 1)));
end