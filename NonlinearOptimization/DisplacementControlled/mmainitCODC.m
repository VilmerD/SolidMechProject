nelm = numel(x0);
xk = x0;
xkm1 = x0;
xkm2 = x0;
xmin = 1e-3*ones(nelm , 1);
xmax = 1e+0*ones(nelm , 1);
low = xmin;
upp = xmax;

maxoutit = 1000;
minoutit = 400;
m = 1;              % Number of constraints
n = nelm;           % Number of dimentions
epsimin = 1e-6;     
a0 = 1;
a = 0;              % Size should match m, and be zero
c = 1e3;            % To make sure y = 0 in any optimal solution
d = 1;              % See above

alph = 1;
alphtol = 3e-7;
