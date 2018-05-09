function node = fun_covProd(n1, n2)
%FUN_COVPROD A product of two covariances.

  node = Node('covProd', [n1.hypInit n2.hypInit], [n1.hypBounds n2.hypBounds]);
  node.setChildren({n1, n2});

end

