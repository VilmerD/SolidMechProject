classdef LinearSolver < handle
    properties
        % Storing old factorizations and solutions
        Rold;
        Aold;
        
        % CA
        its_since_fact;
        force_factorization = 1;
        maxits;
        s;
        
        %
        approx_sens = 0;
        
        %
        statistics;
    end
    
    methods
        function obj = LinearSolver(maxits, s)
            obj.maxits = maxits;
            obj.its_since_fact = maxits;
            
            obj.s = s;
            obj.statistics = struct('ncalls',         0, ...
                                    'factorizations', [], ...
                                    'z',              [], ...
                                    'residuals',      [], ...
                                    'designs',        [], ...
                                    'designUpdate',   []);
        end
        
        % Solves the equilibrium equations given by K, f, and bc
        function [f, x] = solveq(obj, K, f, bc)
            if nargin < 4
                x = K\f;
                f = zeros(length(x), 1);
            else
                obj.statistics.ncalls = obj.statistics.ncalls + 1;
                % Initializing the dofs that correspond to known displacements
                % (bck) and known forces(ddk)
                bck = bc(:, 1);
                xp = bc(:, 2);

                ndof = length(f);
                ddk = 1:ndof;
                ddk(bck) = [];

                % Initializes the stiffness matrix and extracts the submatrices
                % corresponding to solving the system
                % [A B][u] = [g]    u unknown,  g known
                % [C D][v] = [h]    v known,    h unknown
                [A, B, C, D] = extractSubmatrices(K, bck, ddk);

                % Solve the reduced system Au = b
                g = f(ddk);
                v = xp;
                b = g - B*v;

                % If direct solvers is used maxits is 1, otherwise CA is used
                if obj.force_factorization || obj.its_since_fact >= obj.maxits
                    fprintf('\n(LS) Factorizing matrix.')
                    obj.its_since_fact = 1;
                    obj.statistics.factorizations(obj.statistics.ncalls) = 1;
                    obj.Aold = A;

                    % Try to do a cholesky, if the matrix is neg def do lu
                    % instead. 
                    try 
                        R = chol(A);
                        obj.Rold = R;
                        u = R\(R'\b);
                    catch
                        fprintf('K < 0: Computing LU.')
                        u = A\b;
                    end

                % Reusing previous factorization R of the submatrix A
                else
                    fprintf('\n(LS) Reusing factorization.')
                    R = obj.Rold;
                    dA = A - obj.Aold;

                    % Generating basis vectors
                    V = CA(R, A, dA, b, obj.s);

                    % Projecting b onto the basis for the solution
                    z = V'*b;
                    u = V*z;

                    obj.statistics.z(:, obj.statistics.ncalls) = z;
                end
                % Compuing the reaction forces
                h = C*u + D*v;
                f(bck) = h;

                % Assembling final solution vector
                x = zeros(ndof, 1);
                x(bck) = v;
                x(ddk) = u;
            end
        end
        
        function l = approxSens(A, dA, R, f, Rb, y)
            l = zeros(length(f), 1);
            df = f - A*Rb*y;
            l(:, end) = R\(R'\(2*y(end)*df));
            
            for j = obj.s-1:-1:1
                fp = 2*y(j)*df - dA*l(:, j + 1);
                l(:, j) = R\(R'\fp);
            end
        end
        
        function recordUpdate(obj, z)
            st = obj.statistics;
            st.designUpdate = [st.designUpdate st.ncalls];
            st.designs = [st.designs z];
            st.residuals = [st.residuals, 0];
            obj.statistics = st;
            
            obj.its_since_fact = obj.its_since_fact + 1;
            obj.force_factorization = 0;
        end
    end
end