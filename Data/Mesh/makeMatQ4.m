function [ex, ey, coord, dof, edof, enod] = makeMatQ4(res)
nx = res(1);
ny = res(2);
nelm = nx*ny;
nnod = (nx + 1)*(ny + 1);
n = 1:nnod;

edge1 = 1           :1          :(nx + 1);
edge2 = (nx + 1)    :(nx + 1)   :nnod;
edge3 = nnod - nx   :1          :nnod;
edge4 = 1           :(nx + 1)   :nnod - nx;
n1 = setdiff(n, union(edge2, edge3))';
n2 = setdiff(n, union(edge3, edge4))';
n3 = setdiff(n, union(edge4, edge1))';
n4 = setdiff(n, union(edge1, edge2))';
enod = [(1:nelm)' n1 n2 n3 n4];

dof = [1:2:2*nnod; 2:2:2*nnod]';
dofx = dof(:, 1);   edofx = dofx(enod(:, 2:end));
dofy = dof(:, 2);   edofy = dofy(enod(:, 2:end));
edof = (1:nelm)';
edof(:, 2:2:9) = edofx;
edof(:, 3:2:9) = edofy;

xcoord = repmat(linspace(0, 1, (nx + 1)), 1, (ny + 1));
ycoord = reshape(repmat(linspace(0, 1, (ny + 1)), (nx + 1), 1), 1, []);
coord = [xcoord', ycoord'];
ex = xcoord(enod(:, 2:end));
ey = ycoord(enod(:, 2:end));
end