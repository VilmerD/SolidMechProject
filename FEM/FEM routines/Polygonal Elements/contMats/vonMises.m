function seff = vonMises(s, mpara)
sxx = s(1); syy = s(2); sxy = s(3);
szz = mpara(2)*(sxx + syy);

seff = sqrt(sxx ^ 2 + syy ^ 2 + szz ^ 2 - sxx*syy - sxx*szz - syy*szz ...
    + 3*sxy^2);
end