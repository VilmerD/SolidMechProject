function lambda = stretch1D(ec, ed)
% Computes the stretch in 1D
x0 = ec(:, 1) - ec(:, 2);
dx = ed(1:3) - ed(4:6);
dx = x0 + dx;

lambda = norm(x0 + a) / norm(x0);
end