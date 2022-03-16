function [u, P, ef, es, NSTEP] = NRDCoptstep(K, r, sfun, bc, u0, nmax)
% Wraps NRDC
solution_found = false;
u0 = [zeros(size(u0, 1), 1) u0];
RESTARTS_MAX = size(u0, 2);
nrestarts = 0;

% Restarts NRDC at steps u0
while ~solution_found
    % Get starting point for iteration
    uold = u0(:, RESTARTS_MAX - nrestarts);
    
    % Attempt a solution
    [u, P, ef, es, flag, NSTEP, N_INNER, relres] = ...
        NRDC(K, r, sfun, bc, uold, nmax); 
    
    % Check results
    if flag == 0
        % Solution found! Get out!
        solution_found = true;
    
    elseif flag == 1 && nrestarts < RESTARTS_MAX - 1
        % Restart with other initial guess
        nrestarts = nrestarts + 1;
        
    elseif nrestarts == RESTARTS_MAX
        % Try with a smaller dz
        
        
    end
    
end

% If solution is still not found there is a problem...
if ~solution_found
    errorStruct.message = 'Couldnt solve problem';
    error(errorStruct);
end

% Log results
printform = '%5i %5i %5i %5.0e';
fprintf(printform, NSTEP, N_INNER, nrestarts, relres);
end