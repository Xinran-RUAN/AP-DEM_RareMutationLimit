function [dm, dp] = one_sided_derivatives(u, dtheta, method)
%ONE_SIDED_DERIVATIVES D^-u and D^+u, or WENO5 one-sided reconstructions.
if nargin < 3 || isempty(method)
    method = 'first-order';
end
u = u(:).';
switch lower(method)
    case {'first-order','godunov','lf'}
        dm = (u - [u(end), u(1:end-1)]) / dtheta;
        dp = ([u(2:end), u(1)] - u) / dtheta;
    case 'weno5'
        dm = src.weno5_left_derivative(u, dtheta);
        dp = src.weno5_right_derivative(u, dtheta);
    otherwise
        error('Unknown derivative method "%s".', method);
end
end
