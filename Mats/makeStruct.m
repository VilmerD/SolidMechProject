load(geomfile)
vert_left = find(abs(coord(:, 1) - 0) < 1e-6);
dof_left = dof(vert_left, :);   dof_left = dof_left(:);
vert_right = find(abs(coord(:, 1) - 0.6) < 1e-6);
dof_right = dof(vert_right, :); dof_right = dof_right(:);
dofs_zerodisp = [dof_left; dof_right];
z = zeros(size(dofs_zerodisp));

% dx = 6e-3;  dy = 3e-3;
dx = -1;     dy = -1;
is_horzcenter = (abs(coord(:, 1) - 0.3) < dx);
is_vertcenter = (abs(coord(:, 2) - 0.1) < dy);
verts_disp = find(is_horzcenter.*is_vertcenter);
dofs_disp = dof(verts_disp, 2);
nz = ones(size(dofs_disp));

bc = [dofs_zerodisp, z
      dofs_disp,    nz];

[nelm, ~] = size(edof);
ndof = max(edof(:));
F = zeros(ndof, 1);
F(verts_disp) = 1;
xcoord = coord(:, 1);   ex = xcoord(enod(:, 2:end));
ycoord = coord(:, 2);   ey = ycoord(enod(:, 2:end));
newfile = [geomfile(1:end-4), 'New.mat'];
save(newfile, 'F', 'bc', 'coord', 'dof', 'edof', 'enod', 'nelm', ...
    'ndof', 'ex', 'ey')