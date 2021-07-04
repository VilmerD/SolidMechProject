function seff = vonMises(s, mpara)
if iscell(s)
    seff = zeros(numel(s), 1);
    for i = 1:numel(s)
        si = s{i};
        sxx = si(1); syy = si(2); sxy = si(3);
        szz = mpara(2)*(sxx + syy);
        
        seff(i) = sqrt(sxx ^ 2 + syy ^ 2 + szz ^ 2 - sxx*syy - ...
            sxx*szz - syy*szz + 3*sxy^2);
    end
    seff = mean(seff);
else
    sxx = s(1); syy = s(2); sxy = s(3);
    szz = mpara(2)*(sxx + syy);
    
    seff = sqrt(sxx ^ 2 + syy ^ 2 + szz ^ 2 - sxx*syy - ...
        sxx*szz - syy*szz + 3*sxy^2);
end
end