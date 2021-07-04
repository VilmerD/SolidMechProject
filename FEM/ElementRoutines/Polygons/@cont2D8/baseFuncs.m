function dN = baseFuncs()
% Derivatives of base functions
dndxi = @(x, e) [-1/4  *(-1 + e)    *(2*x + e)
                 1     *( 0 + x)    *( -1 + e)
                 1/4   *(-1 + e)    *(  e - 2*x)
                 -1/2  *( 1 + e)    *( -1 + e)
                 1/4   *( 1 + e)    *(2*x + e)
                 -1    *( 0 + x)    *(  1 + e)
                 -1/4  *( 1 + e)    *(  e - 2*x)
                 1/2   *( 1 + e)    *( -1 + e)]';
          
dndeta = @(x, e) [-1/4  *(-1 + x)   *(  x + 2*e)
                  1/2   *( 1 + x)   *( -1 + x)
                  1/4   *( 1 + x)   *(2*e - x)
                  -1    *( 0 + e)   *(  1 + x)
                  1/4   *( 1 + x)   *(  x + 2*e)
                  -1/2  *( 1 + x)   *( -1 + x)
                  -1/4  *(-1 + x)   *(2*e - x)
                  1     *( 0 + e)   *( -1 + x)]';

% Integrating derivatives
x = [-0.774596669241483
     0
     0.774596669241483];
numint = 3;
npoints = 8;

dN = zeros(2*numint^2, npoints);
for k = 1:numint
    for i = 1:numint
        ind = (k - 1)*numint + i;
        dN(ind, :)              = dndxi(x(i), x(k));
        dN(ind + numint^2, :)   = dndeta(x(i), x(k));
    end
end
end