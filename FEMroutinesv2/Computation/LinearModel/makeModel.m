function model = makeModel(geomfile)
load(geomfile, 'edof', 'coord', 'enod', 'F', 'bc');

E = 210e9;
nu = 0.3;
rho = 7800;
t = 1e-3;

xcoord = coord(:, 1);   ex = xcoord(enod(:, 2:end));
ycoord = coord(:, 2);   ey = ycoord(enod(:, 2:end));
model = contModel(edof, ex, ey, [E, nu], t, rho, bc, F);
end