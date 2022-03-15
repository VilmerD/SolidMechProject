function dw = sCAEEONmod(B, psi, w, K, M, Kold, R, psi0, edof, np, K0, M0)
% sCAEEON computes the sensitivity of the eigenvalues when
% eigenmode orthogonalized CA is used to reduce the eigenproblem

ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];
Kff = K(nf, nf);
Mff = M(nf, nf);

dKff = Kff - Kold(nf, nf);

% Preallocating memory
ne = numel(w);
Zs = cell(ne, 1);
Us = cell(ne, 1);
Ts = cell(ne, 1);

QUs = cell(ne, 1);
WPs = cell(ne, 1);
uns = cell(ne, 1);

mpsi0 = M*psi0;
% Computing adjoint vectors
for k = 1:ne
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
        wspj = ws'*psi0(nf, j);
        WP(sk, j) = wspj;
        qs = qs - mpsi0(nf, j)*wspj;
    end
    
    us = Uk(nf, sk);
    qsus = qs'*us;
    
    % Compute M*u in advance to save some time
    mu = M*Uk;
    mus = mu(nf, sk);
    usn = sqrt(us'*mus);
    Zk = zeros(ndof, sk);
    Zk(nf, sk) = R\(R'\((usn^2*qs - qsus*mus)))/usn^3;
    
    QU = zeros(sk, 1);
    QU(sk) = qsus;
    un = zeros(sk, 1);
    un(sk) = usn;
    
    % Stepping backward
    for i = (sk-1):-1:1
        wi = -2*yk(i)*dff;
        qi = wi - dKff*Zk(nf, i+1);
        
        for j = 1:(k - 1)
            wipj = wi'*psi0(nf, j);
            qi = qi - mpsi0(nf, j)*wipj;
            WP(i, j) = wipj;
        end
        
        ui = Uk(nf, i);
        qiui = qi'*ui;
        uin = sqrt(ui'*mu(nf, i));
        Zk(nf, i) = R\(R'\((uin^2*qi - (mu(nf, i))*(qiui))))/uin^3;
        
        QU(i) = qiui;
        un(i) = uin;
    end
    Zs{k} = Zk;
    Us{k} = Uk;
    Ts{k} = Tk;
    QUs{k} = QU;
    uns{k} = un;
    WPs{k} = WP;
end

% Computing sensitivity on the element level
nelm = size(edof, 1);
dl = zeros(nelm, 1);
for elm = 1:nelm
    % Preparing matrices
    dofs = edof(elm, 2:end);
    psi0e = psi0(dofs, :);
    psie = psi(dofs, :);
    
    ke = K0{elm};
    me = M0{elm};
    
    % Computing sensitivity
    for k = 1:ne
        Zelm = Zs{k}(dofs, :);
        Uelm = Us{k}(dofs, :);
        Telm = Ts{k}(dofs, :);
        QU = QUs{k};
        WP = WPs{k};
        un = uns{k};
        
        psiek = psie(:, k);
        psi0ek = psi0e(:, k);
        dlelm = psiek'*(ke - w(k)^2*me)*psiek;
        
        zterm = - (Zelm(:, 1))'*me*psi0ek;
        
        u1elm = Uelm(:, 1);
        qterm = 1/2*QU(1)*(u1elm'*me*u1elm)/un(1)^3;
        for i = 2:sk
            zterm = zterm + (Zelm(:, i))'*ke*Telm(:, i-1);
            
            uielm = Uelm(:, i);
            qterm = qterm + 1/2*QU(i)*(uielm'*me*uielm)/un(i)^3;
        end
        
        if k > 1
            wterm = sum(WP.*(Telm'*me*psi0e(:, 1:(k - 1))), 'all');
        else
            wterm = 0;
        end
        
        dl(elm, k) = dlelm + zterm + qterm + wterm;
    end
end
dw = dl./(2*w');
end