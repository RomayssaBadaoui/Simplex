
function col = MRTD(a, b)

% The Maximum Ratio Test performed on vectors a and b.
% This function is called from within the function dsimplex.
% Output parameter:
% col - index of the pivot column. 

m = length(a);
c = 1:m;
a = a(2:m); %(:)
b = b(2:m); %(:)
l = c(b < 0);
[mi, col] = max(a(l)./b(l));
col = l(col);
