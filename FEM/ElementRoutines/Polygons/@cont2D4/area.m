function A = area(dJ)
w = [1
     1];

ind = [1 1
       2 1
       2 2
       1 2];
A = 0;
for k = 1:4
    A = A + dJ{k} * w(ind(k, 1))*w(ind(k, 2));
end
end