function F = animate2D(Ex, Ey, Ed, pos, fext, tmax, magn, ax)
[~, n] = size(fext);

dt = tmax/n;
t = 0;
for k = 1:n
    fkext = fext(:, k)/1000;
    edk = Ed{k};    
    eldisp2(Ex, Ey, edk, [1 4 1], magn);
    axis(ax)
    
    t = t + dt;
    format = 'magn: %i \nt: %.2fms \ntot load: %2.fkN';
    tex = sprintf(format, magn, t*1000, fkext);
    H = text(pos(1), pos(2), tex);
    
    drawnow;
    F(k) = getframe(gcf);
    delete(H)
    cla;
end
end