% Modified version of code downloaded from 
%Sagi Katz (2025). Hidden Point Removal (https://www.mathworks.com/matlabcentral/fileexchange/16581-hidden-point-removal), MATLAB Central File Exchange. Retrieved August 9, 2025.

function visiblePtInds=HPR(p,C,param)
% HPR (Hidden Point Removal) Method
% "Direct Visibility of Point Sets", Katz S., Tal
% A. and Basri R., SIGGRAPH 2007, ACM Transactions on Graphics, Volume 26, Issue 3, August 2007.

% Edit: 11/13/2024

[~, dim] = size(p);
p = p - C;
normp = sqrt(sum(p.^2, 2));
R = max(normp)*(10^param);
P = p + 2*(R-normp).*(p./normp);
visiblePtInds = unique(convhulln([P; zeros(1,dim)]));
visiblePtInds(end) = [];
