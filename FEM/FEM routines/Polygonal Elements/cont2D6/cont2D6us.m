function F = cont2D6us(ec, ecinc, Fold)
[dN, ~] = cont2D6N(ec);

function f = Fi(p1, p2, p3, fold)
    defm = zeros(4, 12);
    defm(1:2, 1:2:11) = dN(p1, p2, p3);
    defm(3:4, 2:2:12) = dN(p1, p2, p3);
    finc = defm * ecinc;
    finc = [1 + finc(1) finc(2)     0
            finc(3)     1 + finc(4) 0
            0           0           1];
    f = finc*fold;
    f = [f(1, 1) f(1, 2) f(2, 1) f(2, 2)]';
end

points = [2 1 1
          1 2 1
          1 1 2] / 3;
      
F = cell(3, 1);
for k = 1:3
    fok = Fold{k};
    fok = [fok(1)    fok(2)   0
           fok(3)    fok(4)   0
           0          0       1];
    p1 = points(k, 1);  p2 = points(k, 2);  p3 = points(k, 3);
    F{k} = Fi(p1, p2, p3, fok);
end
end