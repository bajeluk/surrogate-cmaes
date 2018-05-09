classdef Node < handle
  %NODE A node in a covariance expression tree.

  properties
    covFcn
    hypInit
    hypBounds
    children
    nHyp
  end

  methods (Access = public)
    function obj = Node(covFcn, hypInit, hypBounds)
      obj.covFcn = covFcn;
      obj.hypInit = hypInit;
      obj.hypBounds = hypBounds;
      obj.children = {};
      obj.nHyp = length(hypInit);
    end

    function setChildren(obj, children)
      obj.children = children;
    end
  end

end

