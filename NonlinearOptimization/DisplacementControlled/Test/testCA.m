%% Running tests
maxits = [4 10 16];
nbasis = [6 9 12];
tests = [distribute(maxits, nbasis) [0; 0]];
ntests = numel(tests)/2;
for k = 1:ntests
    maxits = tests(1, k);
    nbasis = tests(2, k);
	solver = LinearSolver(maxits, nbasis);
    DC;
end