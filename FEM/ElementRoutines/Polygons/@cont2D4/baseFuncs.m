function dN = baseFuncs()
% Derivatives of base functions
dndxi = @(xi, e) [-1/4*(1 - e)
                   1/4*(1 - e)
                   1/4*(1 + e)
                  -1/4*(1 + e)];
              
dndeta = @(xi, e) [-1/4*(1 - xi)
                  -1/4*(1 + xi)
                   1/4*(1 + xi)
                   1/4*(1 - xi)];
               
% Integrating derivatives
x = [-0.577350269189626
     0.577350269189626];
numint = 2;
npoints = 4;

dN = zeros(2*numint^2, npoints);
for k = 1:numint
    for i = 1:numint
        ind = (k - 1)*numint + i;
        dN(ind, :)              = dndxi(x(i), x(k));
        dN(ind + numint^2, :)  = dndeta(x(i), x(k));
    end
end
end

