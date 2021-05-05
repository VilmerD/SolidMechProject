classdef LinearSolver < handle
    properties
        % Storing old factorizations and solutions
        Rold;
        Kold;
        nsteps;
        
        % CA
        its_since_fact;
        force_factorization = 1;
        maxits;
        s;
        
        %
        cothmax;
        
        %
        statistics;
    end
    
    methods
        function obj = LinearSolver(maxits, s, cothmax)
            obj.maxits = maxits;
            obj.its_since_fact = maxits;

            obj.s = s;
            obj.statistics = struct('ncalls',         0, ...
                                    'factorizations', [], ...
                                    'z',              [], ...
                                    'residuals',      [], ...
                                    'designs',        [], ...
                                    'designUpdate',   []);
            if nargin < 3
                cothmax = 1-1e-3;
            end
            obj.cothmax = cothmax;
        end
        
        % Solves the equilibrium equations given by K, f, and bc
        function [f, x] = solveq(obj, K, f, bc, n)
            obj.statistics.ncalls = obj.statistics.ncalls + 1;

            % Initializing the dofs that correspond to known displacements
            % (bck) and known forces(ddk)
            np = bc(:, 1);
            xp = bc(:, 2);

            ndof = length(f);
            nf = 1:ndof;
            nf(np) = [];

            % Initializes the stiffness matrix and extracts the submatrices
            % corresponding to solving the system
            % [A B][u] = [g]    u unknown,  g known
            % [C D][v] = [h]    v known,    h unknown
            [Kff, Kfp, Kpf, Kpp] = extractSubmatrices(K, np, nf);

            % Solve the reduced system Au = b
            ff = f(nf);
            up = xp;
            b = ff - Kfp*up;

            % If direct solvers is used maxits is 1, otherwise CA is used
            if obj.force_factorization || obj.its_since_fact >= obj.maxits
                fprintf('\n(LS) Factorizing matrix, ')
                obj.its_since_fact = 1;
                obj.statistics.factorizations(obj.statistics.ncalls) = 1;
                obj.Kold{n} = Kff;

                % Try to do a cholesky, if the matrix is neg def do lu
                % instead.
                fprintf('trying Cholesky... ')
                try 
                    R = chol(Kff);
                    fprintf('Success. ')
                    obj.Rold{n} = R;
                    uf = R\(R'\b);
                catch
                    % Does not really work yet as i cannot save both U and
                    % L, so this just makes the method fail if the
                    % factorization is to be reused nest step.
                    
                    % maybe using a struct containing information of the
                    % factorization used could solve it.
                    fprintf('Failed, computing LU. ')
                    [L, U] = lu(Kff);
                    uf = U\(L\b);
                end

            % Reusing previous factorization R of the submatrix A
            else
                fprintf('\n(LS) Reusing factorization.')
                obj.statistics.factorizations(obj.statistics.ncalls) = 0;
                R = obj.Rold{n};
                dK = Kff - obj.Kold{n};

                % Generating basis vectors
                V = CA(R, Kff, dK, b, obj.s);

                % Projecting b onto the basis for the solution
                z = V'*b;
                uf = V*z;

                obj.statistics.z(:, obj.statistics.ncalls) = z;
            end
            % Compuing the reaction forces
            fp = Kpf*uf + Kpp*up;
            f(np) = fp;

            % Assembling final solution vector
            x = zeros(ndof, 1);
            x(np) = up;
            x(nf) = uf;
        end
        
        % Checks if the a factorization should be forced
        function checkAngle(obj, znew)
            if obj.statistics.ncalls > 0
                zold = obj.statistics.designs(:, end);
                coth = zold'*znew/(norm(zold)*norm(znew));
                if coth < obj.cothmax
                    obj.force_factorization = 1;
                end
            else
                obj.force_factorization = 1;
            end
        end
        
        % Records data of the design update
        function recordUpdate(obj, z)
            st = obj.statistics;
            facts = st.factorizations;
            st.designUpdate = [st.designUpdate st.ncalls];
            st.designs = [st.designs z];
            st.residuals = [st.residuals, 0];
            obj.statistics = st;
            
            obj.its_since_fact = obj.its_since_fact + 1;
            obj.force_factorization = 0;
        end
        
    end
end