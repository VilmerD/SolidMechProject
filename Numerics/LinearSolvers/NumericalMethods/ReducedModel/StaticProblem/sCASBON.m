function dc = sCASBON(B, f, K, Kold, R, edof, np, K0)
% sCASBON computes the sensitivities of the compliance when CA is
% considered
nelm = size(edof, 1);
ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];

Kff = K(nf, nf);
dKff = Kff - Kold(nf, nf);
ff = f(nf);

Uf = B{1};
Tf = B{2};
Rf = B{3};
Vf = B{4};

s = size(Uf, 2);
U = zeros(ndof, s);
T = zeros(ndof, s);
Rb = zeros(ndof, s);
V = zeros(ndof, s);
U(nf, :) = Uf;
T(nf, :) = Tf;
Rb(nf, :) = Rf;
V(nf, :) = Vf;

% Computing adjoint vectors
ps = -2*ff*(Vf(:, s)'*ff);

rs = Rb(nf, s);
nrs = sqrt(rs'*Kff*rs);
psrs_nrs = ps'*rs/nrs^3;
ws = (ps/nrs - (Kff*rs)*(psrs_nrs));

PR_NR = zeros(s, 1);
PR_NR(s) = psrs_nrs;

qs = ws;
WV = zeros(s, s);
Kv = zeros(ndof, s);
for j = 1:(s - 1)
    vj = V(nf, j);
    wsvj = ws'*vj;
    Kvj = Kff*vj;
    qs = qs - Kvj*wsvj;
    WV(s, j) = wsvj;
    Kv(nf, j) = Kvj;
end
us = U(nf, s);
nus = sqrt(us'*Kff*us);
qsus_nus = qs'*us/nus^3;
zs = (R\(R'\(qs/nus - (Kff*us)*(qsus_nus))));

QU_NU = zeros(s, 1);
QU_NU(s) = qsus_nus;

P(nf, s) = ps;
W(nf, s) = ws;
Q(nf, s) = qs;
Z(nf, s) = zs;
for i = (s-1):-1:1
    pi = -2*ff*(Vf(:, i)'*ff);
    vi = Vf(:, i);
    for j = (i + 1):s
        wj = W(nf, j);
        wjvi = wj'*vi;
        tj = T(nf, j);
        Ktj = Kff*tj;
        pi = pi - wjvi*Ktj - (vi'*Ktj)*wj;
        WV(j, i) = wjvi;
    end
    P(nf, i) = pi;
    
    ri = Rb(nf, i);
    nri = sqrt(ri'*Kff*ri);
    piri_nri = pi'*ri/nri^3;
    wi = (pi/nri - (Kff*ri)*(piri_nri));
    W(nf, i) = wi;
    PR_NR(i) = piri_nri;
    
    zip1 = Z(nf, i + 1);
    qi = wi - dKff*zip1;
    for j = 1:(i-1)
        vj = Vf(:, j);
        wivj = wi'*vj;
        qi = qi - Kv(nf, j)*wivj;
        WV(i, j) = wivj;
    end
    Q(nf, i) = qi;
    
    ui = U(nf, i);
    nui = sqrt(ui'*Kff*ui);
    qiui_nui = qi'*ui/nui^3;
    QU_NU(i) = qiui_nui;
    if i > 1
        zi = (R\(R'\(qi/nui - (Kff*ui)*qiui_nui)));
        Z(nf, i) = zi;
    end
end

% Index matrix to extract element data
index_matrix = reshape(0:ndof:(s-1)*ndof, 1, 1, s) + ...
    repmat(edof(:, 2:end)', 1, 1, s);
Zs = Z(index_matrix);
Us = U(index_matrix);
Ts = T(index_matrix);
Rs = Rb(index_matrix);
Vs = V(index_matrix);
% Computing sensitivities
dc = zeros(nelm, 1);
for elm = 1:nelm
    ke = K0{elm};
    
    Zelm = Zs(:, elm, :);
    Uelm = Us(:, elm, :);
    Telm = Ts(:, elm, :);
    Relm = Rs(:, elm, :);
    Velm = Vs(:, elm, :);
    QU_NUe = QU_NU;
    PR_NRe = PR_NR;
    WVe = WV;
    
    dcelm = 0;
    for i = 1:s
        tie = Telm(:, i);
        if i < s
            dcelm = dcelm + Zelm(:, i + 1)'*ke*tie;
        end
        
        for j = 1:(i - 1)
            dcelm = dcelm + (tie'*ke*Velm(:, j))*WVe(i, j); 
        end
        
        uie = Uelm(:, i);
        dcelm = dcelm + 1/2*QU_NUe(i)*(uie'*ke*uie);
        
        rie = Relm(:, i);
        dcelm = dcelm + 1/2*PR_NRe(i)*(rie'*ke*rie);
    end
    dc(elm) = dcelm;
end