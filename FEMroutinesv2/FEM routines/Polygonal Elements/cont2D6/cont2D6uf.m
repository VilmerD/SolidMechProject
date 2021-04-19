function fe = cont2D6uf(ec, t, stress)
[dN, J] = cont2D6N(ec);

function Bk = Bmat(p1, p2, p3)
    Bk = zeros(3, 12);
    dNk = dN(p1, p2, p3);
    Bk(1, 1:2:11) = dNk(1, :);
    Bk(2, 2:2:12) = dNk(2, :);
    Bk(3, 2:2:12) = dNk(1, :);
    Bk(3, 1:2:11) = dNk(2, :);
end

Hi = 1/3;
pts = [2 1 1
       1 2 1
       1 1 2]/3;
fe = zeros(12, 1);
for k = 1:3
    p1 = pts(k, 1); p2 = pts(k, 2); p3 = pts(k, 3);
    fe = fe + Bmat(p1, p2, p3)'*stress{k}*det(J(p1, p2, p3))*t*Hi;
end

end