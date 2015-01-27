function rm = frmse(x)
  rm = sqrt(sum(x.^2)) / length(x);
end
