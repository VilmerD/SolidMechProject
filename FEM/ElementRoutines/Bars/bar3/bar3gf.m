function fe = bar3gf(ec, ed, es)
% bar3gf(ec, ed, ep) computes the internal force vector given the
% displacement dx, cross-section area A, stress sg, and initial length l0
x0 = ec(:, 1) - ec(:, 2);
l0 = norm(x0);

dx = ed(1:3) - ed(4:6);
dx = x0 + dx;

fe = es / l0 * [dx; -dx];