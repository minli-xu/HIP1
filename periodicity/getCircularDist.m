function d = getCircularDist (c,N)

c2 = c - N;

c1d = dist(c');
c1c2d = pdist2(c,c2);

cmin = min(c1d,c1c2d);
d = nonzeros(triu(cmin));
