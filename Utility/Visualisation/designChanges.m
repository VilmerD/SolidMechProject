function alph = designChanges(X, d)
if nargin == 1
    d = 1;
end
numberdesigns = size(X, 2);
alph = zeros(1, numberdesigns);
for k = 2:numberdesigns
    xkmd = X(:, max(1, k -d));
    xk = X(:, k);
    alph(k) = dot(xkmd/norm(xkmd), xk/norm(xk));
end
end