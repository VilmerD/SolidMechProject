%% Test the filters' derivative
%% SIMP
xmin = 1e-4;
q = 3;
sf = SIMPFilter(xmin, q);
SIMP = @(x) sf.forward(x);
dSIMP = @(x) sf.backward(x);

xx = linspace(0, 1)';
dx = full(diag(dSIMP(xx)));

h = 1e-9;
xxph = xx + h;
dxnum = (SIMP(xxph) - SIMP(xx))/h;

err = norm(dx - dxnum)
rerr = err/norm(dxnum)

%% Mass
mf = MassFilter();
Mass = @(x) mf.forward(x);
dMass = @(x) mf.backward(x);

xx = linspace(0, 1)';
dx = full(diag(dMass(xx)));

h = 1e-9;
xxph = xx + h;
dxnum = (Mass(xxph) - Mass(xx))/h;

err = norm(dx - dxnum)
rerr = err/norm(dxnum)
%% RAMP
xmin = 1e-4;
q = 8;
RAMP = @(x) RAMPFilter.RAMPForward(x, xmin, q);
dRAMP = @(x) RAMPFilter.RAMPBackward(x, xmin, q);

xx = linspace(0, 1)';
dx = full(diag(dRAMP(xx)));

h = sqrt(eps)*100;
xxph = xx + h;
dxnum = (RAMP(xxph) - RAMP(xx))/h;

err = norm(dx - dxnum)
rerr = err/norm(dxnum)
%% Density filter
model = contModel.makeModel('beamFull90x15.mat');
DF = DensityFilter(model.ex, model.ey, model.areas(), 10e-3);

xx = rand(model.nelm, 1);
dx = DF.backward(xx);

h = sqrt(eps)*100;
indx = [10 90 1300 218 555];
dxnum = zeros(model.nelm, numel(indx));
for i = 1:numel(indx)
    ii = indx(i);
    xxph = xx;
    xxph(ii) = xxph(ii) + h;
    
    dxnum(:, i) = (DF.forward(xxph) - DF.forward(xx))/h; 
end

err = sqrt(sum((dx(:, indx) - dxnum).^2, 1))
rerr = err./sqrt(sum((dxnum).^2, 1))

%% Thresholding
b = 100;
eta = 0.5;
H = @(x) HeavisideProjection.H(x, b, eta);
dH = @(x) HeavisideProjection.dH(x, b, eta);

xx = linspace(0, 1)';
dx = full(diag(dH(xx)));

h = sqrt(eps)*100;
xxph = xx + h;
dxnum = (H(xxph) - H(xx))/h;

err = norm(dx - dxnum)
rerr = err/norm(dxnum)

%% Density then SIMP
model = contModel.makeModel('beamFull90x15.mat');
DFilter = DensityFilter(model.ex, model.ey, model.areas(), 10e-3);
DF = @(x) DFilter.forward(x);
dDF = @(x) DFilter.backward(x);

xmin = 1e-9;
q = 3;
SIMP = @(x) SIMPFilter.SIMPForward(x, xmin, q);
dSIMP = @(x) SIMPFilter.SIMPBackward(x, xmin, q);

xx = rand(model.nelm, 1);
dx = dSIMP(DF(xx))*dDF(xx);

h = sqrt(eps)*100;
indx = [10 90 1300 218 555];
dxnum = zeros(model.nelm, numel(indx));
for i = 1:numel(indx)
    ii = indx(i);
    xxph = xx;
    xxph(ii) = xxph(ii) + h;
    
    dxnum(:, i) = (SIMP(DF(xxph)) - SIMP(DF(xx)))/h; 
end

err = sqrt(sum((dx(:, indx) - dxnum).^2, 1))
rerr = err./sqrt(sum((dxnum).^2, 1))

%% SIMP then threshold
xmin = 1e-4;
q = 3;
SIMP = @(x) SIMPFilter.SIMPForward(x, xmin, q);
dSIMP = @(x) SIMPFilter.SIMPBackward(x, xmin, q);

b = 10;
eta = 0.8;
H = @(x) HeavisideProjection.H(x, b, eta);
dH = @(x) HeavisideProjection.dH(x, b, eta);

xx = linspace(0, 1)';
dx = full(diag(dH(SIMP(xx))*dSIMP(xx)));

h = sqrt(eps)*100;
xxph = xx + h;
dxnum = (H(SIMP(xxph)) - H(SIMP(xx)))/h;

err = norm(dx - dxnum)
rerr = err/norm(dxnum)

%% Density then Threshold
model = contModel.makeModel('beamFull90x15.mat');
DFilter = DensityFilter(model.ex, model.ey, model.areas(), 10e-3);
DF = @(x) DFilter.forward(x);
dDF = @(x) DFilter.backward(x);

b = 40;
eta = 0.5;
H = @(x) HeavisideProjection.H(x, b, eta);
dH = @(x) HeavisideProjection.dH(x, b, eta);

xx = rand(model.nelm, 1);
dx = dH(DF(xx))*dDF(xx);

h = sqrt(eps)*100;
indx = [10 90 1300 218 555];
dxnum = zeros(model.nelm, numel(indx));
for i = 1:numel(indx)
    ii = indx(i);
    xxph = xx;
    xxph(ii) = xxph(ii) + h;
    
    dxnum(:, i) = (H(DF(xxph)) - H(DF(xx)))/h; 
end

err = sqrt(sum((dx(:, indx) - dxnum).^2, 1))
rerr = err./sqrt(sum((dxnum).^2, 1))