M = toeplitz([-1 2 zeros(1, 10000)]);
M = M*M' + M'*M;