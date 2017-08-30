X = [0 0; 1 0; 0 1; 1 1]
y = [2; 3; 3; 4]
X1 = [X ones(size(X, 1), 1)]
w = [1 1 2]
X1 * w'
A = [X y ones(size(X, 1), 1)]
l = [1 1 -1 2]'
A * l
M = A' * A
[U,L,W] = eig(M)
L = diag(L)
l = U(:,1);% * Lambda(1,1)
A * l
J = zeros(size(U));
for k = 2:size(U,2)
  lambda = L(k);
  u = U(:,k);% * lambda;
  J = J - u * u' / lambda;
end
S = zeros(size(U));
for i = 1:size(A, 1)
  a = A(i,:);
  Delta = a' * a;
  Delta(end,:) = 0;
  Delta(:,end) = 0;
  S = S + a' * a * (l' * Delta * l);
end
Delta_l = J * S * J;
ld = [l(1:end-2,:); l(end,:)] / l(end-1,:)
f_grad = -eye(size(ld,1), size(l,1)) / l(end-1,:)
f_grad(:, end) = f_grad(:, end-1)
f_grad(:, end-1) = -ld / l(end-1,:)
Delta_ld = f_grad * Delta_l * f_grad'
x = [0 0]
xh = [x 1]
xh * ld
xh * Delta_ld * xh'
ld = -ld
normpdf(2, xh * ld, xh * Delta_ld * xh')