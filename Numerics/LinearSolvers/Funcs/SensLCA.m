function l = SensLCA(A, dA, R, f, Rb, y)
l = zeros(length(f), 1);
df = f - A*Rb*y;
l(:, end) = R\(R'\(2*y(end)*df));

for j = obj.s-1:-1:1
    fp = 2*y(j)*df - dA*l(:, j + 1);
    l(:, j) = R\(R'\fp);
end
end
