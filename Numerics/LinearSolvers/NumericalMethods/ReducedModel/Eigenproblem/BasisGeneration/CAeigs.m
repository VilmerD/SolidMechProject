function [V, D, B] = CAeigs(K, M, n, R, dK, psi0, s, feig)
% CAeigs solves the generalized eigenproblem using a reduced order model
%
% INPUT:
%       K:      Stiffness matrix                                    (n x n)
%       M:      Mass matrix                                         (n x n)

% OUTPUT:
%       V:      Eigenvectors                                        (n x k)
%       D:      Eigenvalues                                         (k x k)

nelm = size(K, 1);
V = zeros(nelm, n);
D = zeros(n, n);
B = cell(n, 2);
for k = 1:n
    % Choose basis generation method
    if strcmp(feig, 'cae')
        Bk = CAE(R, dK, K, M, psi0(:, k), s);
    elseif strcmp(feig, 'caeeon') || strcmp(feig, 'caeeonlazy')
        Bk = CAEEON(R, dK, K, M, [psi0(:, k) V], s);
    elseif strcmp(feig, 'caeeonmod')
        Bk = CAEEON(R, dK, K, M, [psi0(:, k) psi0(:, 1:(k - 1))], s);
    elseif strcmp(feig, 'caers')
        Bk = CAERS(R, dK, K, M, psi0(:, k), s);
    else
        errorStruct.message = sprintf('CA-method not available: %s', feig);
        error(errorStruct);
    end
    
    Vk = Bk{end};
    
    % Compute reduced model
    KR = Vk'*K*Vk;
    MR = Vk'*M*Vk;
    
    [yk, dk] = eigs(KR, MR, 1, 'smallestabs');  % Solve reduced model
    
    yk = yk/sqrt(yk'*MR*yk);                    % Normalize wrt mass matrix
    psik = Vk*yk;                               % Full solution
    
    V(:, k) = psik;                             % Insert solution
    D(k, k) = dk;
    B{k, 1} = Bk;
    B{k, 2} = yk;
end