n = 1000;
r = 5;
X = rand(n, 2) * 2*r - r;

g = -1:2/sqrt(n):1;
[G1,G2] = meshgrid(g);

for r = [5, 10, 100]
  X1 = G1 * r;
  X2 = G2 * r;
  X = [X1(:), X2(:)];
  for i = 1:24
    y = benchmarks(X', i)';

    figure;
    scatter3(X(:,1), X(:,2), y);
    filename = sprintf('figures/f%02d_r%03d.png', i, r);
    print(filename, '-dpng');
    
    figure;
    surf(X1, X2, reshape(y, size(X1)));
    filename = sprintf('figures/f%02d_r%03d_surf.png', i, r);
    print(filename, '-dpng');

    close all;
    i
  end
  r
end