function [w, wmin] = eigenvalues(model, X, ne, p)
nits = size(X, 2);
model.thresholding.update(0);

w = zeros(ne, nits);

for k = 1:nits
    xk = X(:, k);
    model.thresholding.update(k);
    
    xf = model.density_filter.forward(xk);
    K = model.stiffness(xf);
    M = model.mass(xf);
    
    [~, L] = meigenSM(K, M, model.bc(:, 1), ne);
    
    w(:, k) = sqrt(L);
    k
end

wmin = sum(w.^(-p), 1).^(-1/p);
end