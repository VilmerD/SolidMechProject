function Rb = CABasis(R, dA, b, s)
% Computes the CA basis, without normalizing
% First basis vector
u0 = R\(R'\b);
r1 = norm(u0);
u0 = u0/r1;

Rb = [u0];
            
% Finding the other basis vectors
uk = u0;
for k = 2:s
    uk = -R\(R'\(dA*uk));
    Rb(:, k) = uk/r1;
end
end