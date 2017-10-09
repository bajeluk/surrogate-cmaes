for d = [2, 5]
  X = rand(250*d, d);
  y = sum(X, 2);
  
  tic;
  m = fitensemble(X, y, 'Bag', 1, 'Tree', 'type', 'regression');
  yPred = m.predict(X);
  immse(y, yPred)
  toc
  
  tic;
  options = struct;
  options.nTrees = 10;
  m = RandomForestModel(options);
  m = m.trainModel(X, y);
  yPred = m.modelPredict(X);
  immse(y, yPred)
  toc
end