function Fit = evaluateModel(x_test, npoints, model)

if (model.modelType == 1)
%    x_test = (x_test - repmat(model.xmin,1,npoints))./ repmat(model.xmax-model.xmin,1,npoints);
 %   x_test = x_test / model.sigma;
    Fit = -RankSVMFunc(model.xtrainEncoded, x_test, npoints, model.N, ... 
            model.nTrain, model.alphas, model.TwoSigmaPow2,model.invsqrtC, model.Xmean_model, ... 
            model.doEncoding, model.kernel, model.kernelParam1, model.kernelParam2, model.normalize, model.xmin, model.xmax); 

end;

if (model.modelType == 2)
    Fit = SVRFunc(model.xtrainEncoded, x_test, npoints, model.N, ... 
          	model.nTrain, model.alphas, model.TwoSigmaPow2, model.invsqrtC, model.Xmean_model, model.doEncoding, model.bvalue); 
end;

if (model.modelType == 3) 
    Fit = RankStructSVMFunc(model.xtrainEncoded, x_test, npoints, model.N, ... 
          	model.nTrain, model.nConstr, model.alphas, model.TwoSigmaPow2, model.invsqrtC, model.Xmean_model, model.doEncoding);
end;