function si = designAngle(x1, x2)
%DESIGNANGLE computes the angle between designs x1 and x2
if norm(x1 - x2) == 0
    si = 0;
    return
end
x1n = x1/sqrt(x1'*x1);
x2n = x2/sqrt(x2'*x2);

co = dot(x1n, x2n);
si = sqrt(1 - min(co, 1).^2);
end

