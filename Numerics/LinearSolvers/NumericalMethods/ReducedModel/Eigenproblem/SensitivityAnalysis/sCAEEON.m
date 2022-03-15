function dw = sCAEEON(B, psi, w, K, M, Kold, R, psi0, edof, np, K0, M0)
% sCAEEON computes the sensitivity of the eigenvalues when
% eigenmode orthogonalized CA is used to reduce the eigenproblem
ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];
Kff = K(nf, nf);
Mff = M(nf, nf);

dKff = Kff - Kold(nf, nf);

% Preallocating memory
nelm = size(edof, 1);
npts = size(edof, 2) - 1;
edof_red = edof(:, 2:end)';
ne = numel(w);
sguess = size(B{1, 2}, 1);
Zs = zeros(npts, nelm, ne, sguess);
Us = Zs;
Ts = Zs;

QU_UNs = zeros(ne, sguess);
WPs = zeros(ne, sguess, ne);

mpsi = M*psi;
% Computing adjoint vectors
for k = 1:ne
    % Extracting data
    Bk = B{k, 1};
    yk = B{k, 2};
    sk = size(yk, 1);
    Uk = zeros(ndof, sk);
    Tk = zeros(ndof, sk);
    Uk(nf, :) = Bk{1};
    Tk(nf, :) = Bk{2};
    
    psikf = psi(nf, k);
    dff = (Kff*psikf - w(k)^2*Mff*psikf);
    
    % Computing the LAST vector first
    ws = -2*dff*yk(sk);
    
    % last q-adjoint
    qs = ws;
    WP = zeros(sk, k-1);
    for j = 1:(k - 1)
        wspj = ws'*psi(nf, j);
        WP(sk, j) = wspj;
        qs = qs - mpsi(nf, j)*wspj;
    end
    
    mu = M*Uk;
    us = Uk(nf, sk);
    mus = mu(nf, sk);
    usn = sqrt(us'*mus);
    qsus_usn = qs'*us/usn^3;
    
    Zk = zeros(ndof, sk);
    Zk(nf, sk) = R\(R'\((qs/usn - qsus_usn*mus)));
    
    QU_UN = zeros(sk, 1);
    QU_UN(sk) = qsus_usn;
    
    % Stepping backward
    for i = (sk-1):-1:1
        wi = -2*yk(i)*dff;
        qi = wi - dKff*Zk(nf, i+1);
        
        for j = 1:(k - 1)
            wipj = wi'*psi(nf, j);
            qi = qi - mpsi(nf, j)*wipj;
            WP(i, j) = wipj;
        end
        
        ui = Uk(nf, i);
        uin = sqrt(ui'*mu(nf, i));
        qiui_uin = qi'*ui/uin^3;
        Zk(nf, i) = R\(R'\((qi/uin - (mu(nf, i))*(qiui_uin))));
        
        QU_UN(i) = qiui_uin;
    end
    
    % Index matrix to extract element data
    index_matrix = reshape(0:ndof:(sk-1)*ndof, 1, 1, sk) + ...
        repmat(edof_red, 1, 1, sk);
    
    Zs(:, :, k, :) = Zk(index_matrix);
    Us(:, :, k, :) = Uk(index_matrix);
    Ts(:, :, k, :) = Tk(index_matrix);
    QU_UNs(k, :) = QU_UN;
    WPs(k, :, 1:k-1) = WP;
end

% Index matrix to extract element data
index_matrix_full = reshape(0:ndof:(ne-1)*ndof, 1, 1, ne) + ...
    repmat(edof_red, 1, 1, ne);
psielm = psi(index_matrix_full);
psi0elm = psi0(index_matrix_full);
dl = zeros(nelm, ne);
% Computing sensitivity on the element level
for elm = 1:nelm
    % Preparing matrices
    psi0e = psi0elm(:, elm, :);
    psie = psielm(:, elm, :);
    QU_UNse = QU_UNs;
    WPse = WPs;
    we = w;
    
    ke = K0{elm};
    me = M0{elm};
    
    Zelm = Zs(:, elm, :, :);
    Uelm = Us(:, elm, :, :);
    Telm = Ts(:, elm, :, :);
    
    % Computing sensitivity
    for k = 1:ne
        Zelmk = Zelm(:, :, k, :);
        Uelmk = Uelm(:, :, k, :);
        Telmk = Telm(:, :, k, :);
        WPk = WPse(k, :, :);
        
        psiek = psie(:, k);
        psi0ek = psi0e(:, k);
        dlelm = psiek'*(ke - we(k)^2*me)*psiek;
        
        zterm = - (Zelmk(:, 1))'*me*psi0ek;
        qterm = 0;
        wterm = 0;
        for i = 1:sk
            tie = Telmk(:, i);
            if i < sk
                zterm = zterm + (Zelmk(:, i + 1))'*ke*tie;
            end
            
            uielm = Uelmk(:, i);
            qterm = qterm + 1/2*QU_UNse(k, i)*(uielm'*me*uielm);
            
            for j = 1:(k - 1)
                wterm = wterm + (tie'*me*psie(:, j))*WPk(1, i, j);
            end
        end
        
        dl(elm, k) = dlelm + zterm + qterm + wterm;
    end
end
dw = dl./(2*w');
end