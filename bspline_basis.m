function [y,x] = bspline_basis(j,dp,t,x)

if nargin < 4
    x = linspace(t(dp), t(end-dp+1), 100);  % allocate points uniformly
end
y = bspline_basis_recurrence(j,dp,t,x);

function y = bspline_basis_recurrence(j,n,t,x)

y = zeros(size(x));
if n > 1
    b = bspline_basis(j,n-1,t,x);
    dn = x - t(j+1);
    dd = t(j+n) - t(j+1);
    if dd ~= 0  % 0/0 is made zero
        y = y + b.*(dn./dd);
    end
    b = bspline_basis(j+1,n-1,t,x);
    dn = t(j+n+1) - x;
    dd = t(j+n+1) - t(j+1+1);
    if dd ~= 0
        y = y + b.*(dn./dd);
    end
elseif t(j+2) < t(end) 
	y(t(j+1) <= x & x < t(j+2)) = 1;
else
	y(t(j+1) <= x) = 1;
end
