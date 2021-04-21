function [P, u] = uNRDispCont(sys, ndof, bc, n)
zcol = zeros(ndof, 1);
u = zeros(ndof, n);

P = zeros(ndof, n);
P(:, 1) = sys.ru(zcol, zcol);

duk = zeros(ndof, 1);

for k = 2:n+1
    uk = u(:, k-1);
    duk = solveq(sys.Ku(duk), zcol, bc);
    
    u(:, k) = uk + duk;
    P(:, k) = sys.ru(duk, zcol);
    sys.update(duk);
end
