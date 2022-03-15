function savedir = generateFolder(parentdir)
[y, m, d, ~, ~, ~] = datevec(now);
base_folder = sprintf('%s%02i%02i%02i', parentdir, mod(y, 2000), m, d);
numdir = numel(dir([base_folder, '*']));
savedir = sprintf('%s-%i', base_folder, numdir);
mkdir(savedir);
end