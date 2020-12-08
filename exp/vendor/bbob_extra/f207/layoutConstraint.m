function [ issafe ] = layoutConstraint( xSphere, ySphere )
%LAYOUTCONSTRANT Summary of this function goes here
%   Detailed explanation goes here
safe_distance = 1;


issafe = true;

c = combnk(1:length(xSphere),2);

for i=1:length(c(:,1))
    
    d = ((xSphere(c(i,1))-xSphere(c(i,2)))^2+(ySphere(c(i,1))-ySphere(c(i,2)))^2)^0.5;
    
%     fprintf('%f\n', d);
    
    if d<safe_distance
        issafe = false;
        break;
    end
end

end

