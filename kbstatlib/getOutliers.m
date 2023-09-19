function idx = getOutliers(x, thresh)
%% Calculate the outliers in a given dataset.
% Outliers are detected as points lying outside the interquartile range
% (IQR) multiplied by the optional argument "thresh" which defaults to 1.5.
%
% SYNTAX
% idx = getOutliers(x, thresh)
%
% INPUT
%   x           (Nx1 or 1xN double) Input data
%   thresh      (double) Threshold to detect outliers. 
%               
%               OPTIONAL, default = 1.5.
%
% OUTPUT
%   idx         (Nx1 or 1xN logical) Logical array of indices of the 
%               outliers in x.
%
% Author: Kim Joris Boström

if nargin < 2
    thresh = 1.5;
end

Q = quantile(x,  [0.25, 0.75]);
IQR = Q(2) - Q(1);
up =  Q(2) + thresh * IQR;
lo =  Q(1) - thresh * IQR;
idx = ~((x <= up) & (x >= lo));
end