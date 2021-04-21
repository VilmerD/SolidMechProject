function E = strain(defgrad)
F = reshape(defgrad, 2, 2)';    F = [F, [0; 0]; [0 0] 1];
E = (F' * F - eye(3))/2;
E = [E(1, 1) E(2, 2) 2*E(1, 2)]';
end