X = rand(1000, 2);
y = X(:, 1) + randn(size(X, 1), 1) * 0.001;

m = fitlm(X, y, 'constant');
[yPred1, ci1] = m.predict(X);

X1 = generateFeatures(X, 'constant', true);
[~, ~, r, rint] = regress(y, X1);
yPred2 = y - r;
ci2 = yPred2 + rint - r;

modelOptions = struct;
modelOptions.modelSpec = 'constant';
mm = PolynomialModel(modelOptions)
mm = mm.trainModel(X, y);
[yPred3, sd2, ci3] = mm.modelPredict(X);

sum(yPred1 - yPred2)
d = ci1 - ci2