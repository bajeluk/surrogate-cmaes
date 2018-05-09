function node = fun_covSum(n1, n2)
%FUN_COVSUM A sum of two covariances.

  node = Node('covSum', [n1.hypInit n2.hypInit], [n1.hypBounds n2.hypBounds]);
  node.setChildren({n1, n2});

end

