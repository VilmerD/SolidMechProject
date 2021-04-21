function Et = dmat1D(E, lambda)
% Computes the tangential stiffness matrix E
Et = E / 5 * (2 + 3 / lambda^5);