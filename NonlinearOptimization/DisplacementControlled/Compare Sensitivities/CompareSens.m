vq = 0.3;

geomfile = 'beamFullCoarseNew.mat';
load(geomfile);

eltype = '2D4t';
t = 1e-3;
mpara = [1, 0.3];

model = NLCont2D([ex, ey], edof, ndof, mpara, t, eltype, bc, 2);
x0 = ones(nelm, 1)*vq;

fr = 10e-3;
model.fr = fr;
c0 = 32;

maxits = 6;
nb = [2, 4, 12, 16];
lnb = length(nb);

statsvector = [];
for nbasis = nb
    solver = LinearSolver(maxits, nbasis);

    objective = SetupCompareSens(model, solver, vq, x0, c0);

    mmainit;
    mmamain;

    statsvector = [statsvector, solver.statistics];
end
save('stats.mat', 'statsvector')
%%
load 'stats.mat'
figure(1)
for k = 1:lnb
    sk = statsvector(k);
    dck = sk.del_dc;
    dc_an = dck(1:4, :);
    dc_num = dck(5:end, :);
    subplot(1, lnb, k)
    semilogy(abs((abs(dc_an) - abs(dc_num)))')
    xlabel('iteration number')
    ylabel('difference');
    title(['number of basis vectors: ', num2str(nb(k))]);   
    fprintf('number of factorizations: %i\n', ...
        sum(countFactorizations(sk.factorizations, sk.designUpdate)));
end