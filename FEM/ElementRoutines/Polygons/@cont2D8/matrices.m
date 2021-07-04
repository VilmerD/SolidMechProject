function [dJ, H, B] = matrices(ec, dN)
nPoints = 8;
nInt = 9;
dNdxi   = dN(1:nInt, :);
dNdeta  = dN(nInt+1:2*nInt, :);

dJ = cell(nInt, 1);
H = cell(nInt, 1);
B = cell(nInt, 1);
for k = 1:nInt
    dNdxik = dNdxi(k, :);
    dNdetak = dNdeta(k, :);
    J = [dNdxik     *ec'
         dNdetak    *ec'];
    dJ{k} = det(J);
    
    dNcurr = J \ [dNdxik; dNdetak];
    
    Blk = zeros(3, nPoints*2);
    Blk(1, 1:2:nPoints*2) = dNcurr(1, :);
    Blk(2, 2:2:nPoints*2) = dNcurr(2, :);
    Blk(3, 2:2:nPoints*2) = dNcurr(1, :);
    Blk(3, 1:2:nPoints*2) = dNcurr(2, :);
    B{k} = Blk;
    
    Hk = zeros(4, nPoints*2);
    Hk(1:2, 1:2:nPoints*2) = dNcurr;
    Hk(3:4, 2:2:nPoints*2) = dNcurr;
    H{k} = Hk;
end
end