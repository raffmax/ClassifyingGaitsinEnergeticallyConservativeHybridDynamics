function [f,jacobian] = getJacobianFD(fcn,x,stepSize)
%GETJACOBIANFD Finite Differences on given function
%   Detailed explanation goes here

% get f(x)
f  = fcn(x);
% get dimensions
nF = length(f);
nX = length(x);

jacobian = zeros(nF,nX);
h        = stepSize;
for iCol = 1:nX
       iX       = zeros(nX,1);
       iX(iCol) = 1;
       f_h      = fcn(x+h*iX);
       jacobian(:,iCol) = (f_h-f)/h;
end

end

