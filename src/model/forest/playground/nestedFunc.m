f1 = @(x) sqrt(x.^2);
g = @(x) x.^2;
h = @(x) sqrt(x);
f2 = @(x) h(g(x));

tic;
for i = 1:100
  x = rand(1000000,1);
  y = sqrt(x' * x);
end
toc

tic;
for i = 1:100
  x = rand(1000000,1);
  y = sqrt(x.^2);
end
toc

tic;
for i = 1:100
  x = rand(1000000,1);
  y = f1(x);
end
toc

tic;
for i = 1:100
  x = rand(1000000,1);
  y = f2(x);
end
toc
