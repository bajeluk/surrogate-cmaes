X = rand(10000, 100);

h = Handle(X);
tic
for i = 1:1000
  h = h.set(1, 1, rand());
end
toc

h = Value(X);
tic
for i = 1:1000
  h = h.set(1, 1, rand());
end
toc
