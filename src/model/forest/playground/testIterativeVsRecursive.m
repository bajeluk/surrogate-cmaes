m = 10000;
X = (1:m)';

T = zeros(m, 3);

for n = 1:m
  h = MyHandle(X(1:n, :));
  
  tic;
  s1 = sum(h.X);
  t1 = toc;
  
  tic;
  s2 = h.sumIterative();
  t2 = toc;
  
  tic;
  s3 = h.sumRecursive();
  t3 = toc;
  
  T(n, :) = [t1 t2 t3];
  assert(s1 == s2 && s1 == s3);
end

plot(T(:,1));
hold on;
plot(T(:,2));
hold on;
plot(T(:,3));
hold on;
hold off;
legend('sum', 'sumIterative', 'sumRecursive');