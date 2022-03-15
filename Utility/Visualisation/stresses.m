function [es_full, et_full] = stresses(model, x)
% Compute deformations
xf = model.density_filter.forward(x);
K = model.stiffness(xf);
u = msolveq(K, model.F, model.bc);
enod = model.edof(:, 2:end);

ex = model.ex;
ey = model.ey;
ep = model.ep;
D = hooke(1, model.mpara(1), model.mpara(2));
ed = u(enod);

nelm = model.nelm;
es_full = zeros(nelm, 3);
et_full = zeros(nelm, 3);

for elm = 1:nelm
    [eselm, etelm] = planqs(ex(elm, :), ey(elm, :), ep, D, ed(elm, :), [0; 0]);
    es_full(elm, :) = eselm;
    et_full(elm, :) = etelm;
end

end