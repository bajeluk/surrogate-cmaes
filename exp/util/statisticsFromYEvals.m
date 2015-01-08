function [medians q1 q3 evals] = statisticsFromYEvals(y_evals, maxfunevals, column)

  evals = 1:maxfunevals';

  M = zeros(maxfunevals,length(y_evals));
  for i = 1:length(y_evals)
    if (column == 1)
      M(:,i) = smoothYEvals(y_evals{i}(:,[1 2]), maxfunevals);
    else
      M(:,i) = smoothStatEvals(y_evals{i}(:,[column 2]), maxfunevals);
    end
  end
  Q = quantile(M, [0.25 0.5 0.75], 2);
  q1 = Q(:,1);
  medians = Q(:,2);
  q3 = Q(:,3);
end
