function animate(vidObj, edof, ex, ey, u, P, tend, xc1, yc1, xc2, yc2, r)
[~, n] = size(u);
ed = cell(n, 1);
for k = 1:n
    uk = u(:, k);
    ukf = uk(edof(:, 2:end));
    
    ed{k} = ukf;
end
axi = [0 0.35 -0.25 0.08];
pos = [0.005 0.065];
   
f = figure(2);
set(f, 'position', [50, 50, 600, 600]);
dt = tend/n;

open(vidObj);
for k = 1:n
    cla;
    t = dt*k;
    fkext = norm(P(:, k))/1000;
    edk = ed{k};
    
    hold on;
    eldraw2(ex, ey, [1 2 1]);
    hold on;
    
    eldisp2(ex, ey, edk, [1 4 1], 1);
    hold on;
    
    format = 'magn: %i \nt: %.2fms \ntot load: %.3fkN';
    tex = sprintf(format, 1, t*1000, fkext);
    axis(axi)
    H = text(pos(1, 1), pos(1, 2), tex);
    
    circle(xc1, yc1, r);
    circle(xc2, yc2, r);
    
    axis(axi)
    Fk = getframe(gcf);
    writeVideo(vidObj, Fk);
    delete(H)
    fprintf('%4.0f%%\n', floor(k/n*100));
end
close(vidObj);
end