classdef PDEFilter < Filter
    properties
        l0;
        ls;
    end
    
    methods
        function obj = PDEFilter(nx, ny, w, h, l0, ls)
            [K, T] = assemblePDEFilter(nx, ny, w, h, l0, ls);
            obj.T = T;
            R = chol(K, 'lower');
            obj.R = R;
            
            dy = obj.R'\(obj.R\(obj.T));
            obj.dy = dy;
            
            obj.forward = @(z) F(z);
            obj.backward = @(z) dFdz(z);
        end
        
        function y = F(obj, z)
            y = obj.R'\(obj.R\(obj.T*z));
        end
        
        function dy = dFdz(obj, ~)
            dy = obj.dy;
        end
    end
    
    methods (Static)
        function [KpMpMs, T] = assemblePDEFilter(nx, ny, w, h, l0, ls)
            % Discretizing the space
            [ex, ey, coord, ~, edof] = makeMatQ4Conv(nx, ny);
            coord = [coord(:, 1)*w coord(:, 2)*h];
            Ae = w*h/(nx*ny);
            
            % Creating base matrices and allocating memory for the sparse matrix
            Kb = Ae*l0^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
                [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
            Msb = [2  1; 1  2]'/6*ls;
            Tb = [1  1  1  1]'/4;
            nelm = size(edof, 1);
            I = reshape(kron(edof(:, 2:end), ones(4, 1))', [], 1);
            J = reshape(kron(edof(:, 2:end), ones(1, 4))', [], 1);
            X = reshape(Kb(:)*ones(1, nelm), [], 1);
            KpM = sparse(I, J, X);
            
            Jt = reshape(repmat(1:nelm, 4, 1), [], 1);
            It = reshape(edof(:, 2:end)', [], 1);
            Xt = reshape(Tb*ones(1, nelm), [], 1);
            T = sparse(It, Jt, Xt);
            
            % Surface integral is treated specially
            % Finding boundary dofs for the surface integral
            cx = coord(:, 1);   cy = coord(:, 2);
            W = max(ex(:));     H = max(ey(:));
            So = find(cy == 0); Ea = find(cx == W);
            No = find(cy == H); We = find(cx == 0);
            pairs = [So(1:end-1) So(2:end); Ea(1:end-1) Ea(2:end)
                No(1:end-1) No(2:end); We(1:end-1) We(2:end)];
            
            % Assembling matrix
            Is = reshape(kron(pairs, ones(2, 1))', [], 1);
            Js = reshape(kron(pairs, ones(1, 2))', [], 1);
            dl = sqrt(sum((coord(pairs(:, 1), :) - coord(pairs(:, 2), :)).^2, 2));
            Xs = reshape(Msb(:)*dl', [], 1);
            Ms = sparse(Is, Js, Xs);
            
            KpMpMs = KpM + Ms;
        end
    end
end