function f = countFactorizations(factorizations, designUpdate)
updates = length(designUpdate) - 2;
f = zeros(updates, 1);

for k = 1:updates
    f(k) = sum(factorizations(designUpdate(k):designUpdate(k + 1)));
end