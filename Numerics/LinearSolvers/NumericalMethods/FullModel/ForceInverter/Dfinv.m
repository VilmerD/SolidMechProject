function du = Dfinv(K, L, bc, u, edof, K0)
lamb = msolveq(K, -L, bc);
nelm = size(edof, 1);
du = zeros(nelm, 1);

uelm = u(edof(:, 2:end)');
lambelm = lamb(edof(:, 2:end)');
for elm = 1:nelm
    ke = K0(elm);
    du(elm) = lambelm(:, elm)'*ke*uelm(:, elm);
end

end