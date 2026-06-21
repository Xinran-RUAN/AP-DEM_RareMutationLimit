function [n, D, K, grid, W0, u0] = initial_direct_1d(par)
%INITIAL_DIRECT_1D Initial density obtained from the WKB initial data.
[W0, u0, D, K, grid] = model.initial_wkb_1d(par);
E = exp(u0 / par.eps);
n = bsxfun(@times, W0, E);
end
