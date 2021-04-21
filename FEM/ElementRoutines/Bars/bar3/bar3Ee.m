function Ke = bar3Ee(ec, ed, edprev, ep, es)
dx0 = ec(:, 1) - ec(:, 2);  l0 = norm(dx0);
du = ed(1:3) - ed(4:6);
dx = dx0 + du;

dup = edprev(1:3) - edprev(1:3);
dubar = (du + dup)/2;
dxbar = dx0 + dubar;

outer = dxbar*dx';
I3 = eye(3);

Ke = ep(1) * ep(2) /l0^3 * [outer, -outer; -outer, outer] + ...
    es / l0 * [I3, -I3; -I3, I3];
end