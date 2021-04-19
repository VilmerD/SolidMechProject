function V = CA(R, K, dK, b, s)
% First basis vector and scaling factor r1
u0 = R\(R'\b);
r1 = norm(u0);
u0 = u0/r1;

% First normalized vector
v1 = u0/(u0'*K*u0)^(0.5);
V = zeros(length(b), s);
V(:, 1) = v1;
            
% Finding the other basis vectors
ui = u0;
for i = 2:s
    % Computing next basis vector
    ui = -R\(R'\(dK*ui));
    ri = ui/r1;
    
    % Orthonormalizing the vector with respect to K such that 
    % vi'Kvj = d_ij
    vi = ri;
    for j = 1:i-1
        vj = V(:, j);
        vi = vi - (ri'*K*vj)*vj;
    end
    vi = vi/(vi'*K*vi)^(0.5);
    V(:, i) = vi;
end
end