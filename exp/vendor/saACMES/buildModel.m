function model = buildModel(cur_state,coeff, modelType)

if (modelType == 1)
    model = xacmes_buildModel_RANKSVM(cur_state,coeff);
end;

if (modelType == 2)
    model = xacmes_buildModel_SVR(cur_state,coeff);
end;

if (modelType == 3)
    model = xacmes_buildModel_RANKSTRUCTSVM(cur_state,coeff);
end;