function [dJ, H, B] = matrices(ec, dN)
dNdxi = dN(1:4, :);
dNdeta = dN(5:8, :);

dJ = cell(4, 1);
H = cell(4, 1);
B = cell(4, 1);
for k = 1:4
    dNdxik = dNdxi(k, :);
    dNdetak = dNdeta(k, :);
    Jk = [dNdxik * ec'
          dNdetak * ec'];
    dJ{k} = det(Jk);
    
    dNcurr = Jk \ [dNdxik; dNdetak];
    
    Blk = zeros(3, 8);
    Blk(1, 1:2:7) = dNcurr(1, :);
    Blk(2, 2:2:8) = dNcurr(2, :);
    Blk(3, 2:2:8) = dNcurr(1, :);
    Blk(3, 1:2:7) = dNcurr(2, :);
    B{k} = Blk;
    
    Hk = zeros(4, 8);
    Hk(1:2, 1:2:7) = dNcurr;
    Hk(3:4, 2:2:8) = dNcurr;
    H{k} = Hk;
end
end