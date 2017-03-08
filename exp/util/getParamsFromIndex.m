function [bbParams, sgParams, cmParams, nNonBbobValues, totalCombs] = getParamsFromIndex(id, bbAll, sgAll, cmAll)
% GETPARAMSFROMINDEX Generates struct arrays with parameters accrd. to exper. ID#
  bbNValues = cell2mat(structMap(bbAll, @(x) length(x.values)));
  sgNValues = cell2mat(structMap(sgAll, @(x) length(x.values)));
  cmNValues = cell2mat(structMap(cmAll, @(x) length(x.values)));

  nValues = [bbNValues sgNValues cmNValues];
  totalCombs = prod(nValues);

  assert(id <= totalCombs && id > 0, 'getParamsFromIndex(): the specified index is greater than total # combinations.');

  % Generate a vector with indices of the respective values. The 'nValues' vector
  % defines the number of values in respective fields, for example
  %  nValues =  [3  3  1  1  2  3  1  1  5  2];
  % means that there are 3 different values in the first fiel', 3 in the second field,
  % one in the third field etc.
  paramIDs = getParamIndexVector(id, nValues);

  % The number of different parameter-values combinations except BBOB settings
  nNonBbobValues = prod(nValues((length(bbAll)+1):end));

  pid = 0;

  % BBOB parameters
  bbParams = struct();
  for i = 1:length(bbAll)
    bbParams.(bbAll(i).name) = bbAll(i).values{paramIDs(i)};
  end
  pid = pid + length(bbAll);

  % Surrogate parameters
  sgParams = struct();
  for i = 1:length(sgAll)
    this_name   = sgAll(i).name;
    this_values = sgAll(i).values;
    f_names = strsplit(this_name, {'.', '__'});

    % behave differently if 'modelOpts' is under the cursor:
    % replace sub-settings separated by __ into separate struct
    if (strcmpi(f_names, 'modelopts'))
      level1_structs = {};
      modelOpts_struct = this_values{paramIDs(pid+i)};
      modelOpts_f_names = fieldnames(modelOpts_struct);
      for mo_i = 1:length(modelOpts_f_names)
        mo_f_name = modelOpts_f_names{mo_i};
        f_name_splitted = strsplit(mo_f_name, {'.', '__'});
        if (length(f_name_splitted) > 1)
          if (~ismember(f_name_splitted{1}, level1_structs))
            level1_structs{end+1} = f_name_splitted{1};
            if (~isfield(modelOpts_struct, f_name_splitted{1}))
              modelOpts_struct.(f_name_splitted{1}) = struct();
            end
          end
          modelOpts_struct.(f_name_splitted{1}).(f_name_splitted{2}) = ...
              modelOpts_struct.(mo_f_name);
          modelOpts_struct = rmfield(modelOpts_struct, mo_f_name);
        end
      end
      this_values{paramIDs(pid+i)} = modelOpts_struct;
    end

    switch length(f_names)
    case 1
      sgParams.(this_name) = this_values{paramIDs(pid+i)};
    case 2
      sgParams.(f_names{1}).(f_names{2}) = this_values{paramIDs(pid+i)};
    case 3
      sgParams.(f_names{1}).(f_names{2}).(f_names{3}) = this_values{paramIDs(pid+i)};
    otherwise
      % join the names with underscores
      sgParams.(strjoin(C,'_')) = this_values{paramIDs(pid+i)};
    end
  end
  pid = pid + length(sgAll);

  % CMA-ES parameters
  cmParams = struct();
  for i = 1:length(cmAll)
    cmParams.(cmAll(i).name) = cmAll(i).values{paramIDs(pid+i)};
  end
end
