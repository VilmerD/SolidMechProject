function K = bar3ge(ec, ed, ep, es)
% bar3ge(ec, ed, ep, es) computes the stiffness matrix for a bar element, 
% given the reference positions ec, the displacement ed, material constants
% ep, and the element normal force es
dx0 = ec(:, 1) - ec(:, 2);
du = ed(1:3) - ed(4:6);
dx = dx0 + du;

outer = dx * dx';
l0 = norm(dx0);

I3 = eye(3);
K = ep(1) * ep(2) / l0^3 * [outer, -outer; -outer, outer] + ...
    es / l0 * [I3, -I3; -I3, I3];