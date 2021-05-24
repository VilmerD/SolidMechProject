function X = distribute(X1, X2)
%DISTRIBUTE distributes X1 and X2 such that each permutation of X1 and X2
%           are in X
nx1 = length(X1);   X1 = X1(:);
nx2 = length(X2);   X2 = X2(:);
X = [reshape(repmat(X1', nx2, 1), 1, nx1*nx2); repmat(X2', 1, nx1)];
end

