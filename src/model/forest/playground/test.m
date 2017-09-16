n = 100;
X = rand(n, 2);
y = [ones(n/2, 1); -ones(n/2, 1)];
splitGain = MSESplitGain('constant');
splitGain = splitGain.reset(X, y);

n = 1000000;
y = rand(n, 1);
yPred = rand(n, 1);

for i = 1:1000
end

tic
for i = 1:1000
    sd2 = mean((y-yPred).^2);
end
toc

tic
for i = 1:1000
    r = y-yPred;
    sd2 = r' * r / numel(r);
end
toc

tic
for i = 1:1000
    sd2 = immse(y,yPred);
end
toc
