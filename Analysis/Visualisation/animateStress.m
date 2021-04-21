function animateStress(vidObj, f, edof, ex, ey, u, S, P, ...
    tend, p1, p2, r)
% Plots the object and the stress in it during simulation
%
% Inputs: A bunch. Don't bother
[~, n] = size(u);
ed = cell(n, 1);
for k = 1:n
    uk = u(:, k);
    ukf = uk(edof(:, 2:end));
    
    ed{k} = ukf;
end
axi = [0 0.400 -0.250 0.140];
cax = [1e-4 1]*max(S(:));
text_pos = [0.005 0.125];
   
set(f, 'position', [50, 50, 500, 500]);
dt = tend/n;

open(vidObj);
for k = 1:n
    cla;
    
    Sk = S(:, k); Ske = Sk(edof(:, 2:2:7))';
    Fk = image(k*dt, P(:, k), ed{k}, ex, ey, Ske, cax, text_pos, axi, ...
        p1, r, p2, r);
    
    writeVideo(vidObj, Fk);
end
close(vidObj);
end

function Fk = image(t, Pk, edk, ex, ey, Ske, cax, pos, axi, p1, r1, p2, r2)
    fkext = norm(Pk)/1000;
    
    hold on;
    exd = edk(:, 1:2:5);     x = (ex + exd)';
    eyd = edk(:, 2:2:6);     y = (ey + eyd)';
    
    fill(x, y, Ske, 'Linestyle', 'None');
    colormap(jet);
    colorbar('East');
    caxis(cax)
    
    format = 't: %.2fms \ntot load: %.3fkN';
    tex = sprintf(format, t*1000, fkext);
    axis(axi)
    
    H = text(pos(1, 1), pos(1, 2), tex);
    circle(p1, r1, 'k');
    circle(p2, r2, 'k');
    
    axis(axi)
    xlabel('x-direction [m]')
    ylabel('y-direction [m]')
    Fk = getframe(gcf);
    delete(H);
end