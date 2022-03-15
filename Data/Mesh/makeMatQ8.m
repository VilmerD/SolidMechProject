function [ex, ey, coord, dof, edof] = makeMatQ8(res)
nx = res(1);
ny = res(2);
nelm = nx*ny;
nnod = (2*nx + 1)*(2*ny + 1) - nx*ny;
n = 1:nnod;

% Finding corner nodes
corners = [];
for r = 1:(ny + 1)
    cr = ((3*r-3)*nx + 2*r - 1):2:((3*r-1)*nx + 2*r - 1);
    corners = [corners cr];
end
% Edges are easy to fnd
edges = setdiff(n, corners);
r = 1:(ny + 1);
rm = 1:ny;
south = 1:(2*nx + 1);
east  = [(3*r - 1)*nx + 2*r - 1, 3*rm*nx + 2*rm];
north = (3*nx*ny + 2*ny + 1):nnod;
west  = [1, 1 + setdiff(east, nnod)];
cols = reshape([0 3*rm*nx + 2*rm]' + 2*(1:nx), 1, []);
rows = setdiff(edges, cols);

n1 = setdiff(corners, union(east, north))';
n2 = setdiff(cols, union(east, north))';
n3 = setdiff(corners, union(north, west))';
n4 = setdiff(rows, union(north, west))';
n5 = setdiff(corners, union(west, south))';
n6 = setdiff(cols, union(west, south))';
n7 = setdiff(corners, union(south, east))';
n8 = setdiff(rows, union(south, east))';
enod = [(1:nelm)' n1 n2 n3 n4 n5 n6 n7 n8];

dof = [1:2:2*nnod; 2:2:2*nnod]';
dofx = dof(:, 1);   edofx = dofx(enod(:, 2:end));
dofy = dof(:, 2);   edofy = dofy(enod(:, 2:end));
edof = (1:nelm)';
edof(:, 2:2:17) = edofx;
edof(:, 3:2:17) = edofy;

bigrow = linspace(0, 1, (2*nx + 1));
smallrow = linspace(0, 1, (nx + 1));
xcoord = [repmat([bigrow smallrow], 1, ny), bigrow];
base = [zeros(1, 2*nx + 1) ones(1, nx + 1)/(2*ny)];
ycoord = base + (linspace(0, 1, ny + 1))';
ycoord = reshape(ycoord', 1, []);
ycoord = ycoord(1:nnod);

coord = [xcoord', ycoord'];
ex = xcoord(enod(:, 2:end));
ey = ycoord(enod(:, 2:end));
end