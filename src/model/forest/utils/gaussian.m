X = [];
y = [];
o1 = [-25 -25]
r1 = 25
o2 = [25 25]
r2 = 25
for i = -100:100
  for j = -100:100
    X = [X; i j];
    f = 0;
    x = [i j];
    if norm(x-o1) < r1
      f = sqrt(r1^2 - norm(x-o1)^2);
    elseif norm(x-o2) < r2
      f = sqrt(r2^2 - norm(x-o2)^2);
    end
    y = [y; f];
  end
end
figure;
scatter3(X(:,1), X(:,2), y);
figure;
plot3(X(:,1), X(:,2), y);