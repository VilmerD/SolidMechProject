function newfile = appendTimestamp(file)
[~, m, d, h, mi, s] = datevec(now);
form = [file, '%02i%02i_%02i%02i%02i'];
newfile = sprintf(form, m, d, h, mi, round(s));
end