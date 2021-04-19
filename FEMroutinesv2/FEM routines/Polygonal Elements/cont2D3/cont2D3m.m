function Me = cont2D3m(ec, t, rho)
% Element functions are very simple for this element
A = det([1 1 1; ec])/2;
row1 = [2 0 1 0 1 0];
Me = toeplitz(row1) * rho * t * A / 12;
end
