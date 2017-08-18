n = 100;
X = rand(n, 2);
y = [ones(n/2, 1) -ones(n/2, 1)];
splitGain = SSESplitGain('constant');
splitGain = splitGain.reset(X, y);