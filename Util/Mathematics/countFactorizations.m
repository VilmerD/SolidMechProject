function f = countFactorizations(factorizations, designUpdate)
updates = length(designUpdate);
f = zeros(updates, 1);
dup = 1;
for k = 1:updates
    du = designUpdate(k);
    f(k) = sum(factorizations(dup:du));
    dup = du;
end