%PLOTSCMAESTXTLOG Plots variables from S-CMA-ES log in .txt logfile
function fig = plotScmaesTxtLog(exp_id, fun, dim, id)
  t = readScmaesTxtLog(exp_id, fun, dim, id);

  fig = figure();

  semilogy(t.totEvals, t.sigmaPow2, 'LineWidth', 2);
  hold on;
  semilogy(t.totEvals, t.Dopt, 'LineWidth', 2);
  semilogy(t.totEvals, t.rmseReeval)
  semilogy(t.totEvals, t.diagD_ratio);
  for i = 21:size(t, 2)
    semilogy(t.totEvals, table2array(t(:,i)));
  end

  legend('sigma^2', 'Dopt', 'rmseReeval', 'diagDratio');
end
