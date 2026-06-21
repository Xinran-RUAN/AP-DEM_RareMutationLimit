function [Wn, un, c] = normalize_gauge(W, u, eps, mode)
%NORMALIZE_GAUGE 平移 u 并相应缩放 W，使 n=W exp(u/eps) 不变。
%   注意：小 eps 下若 c/eps 很大，W 会被指数放大。新版 split-if 默认
%   使用 incremental gauge，不再每步强制 max(u)=0。
if nargin < 4 || isempty(mode)
    mode = 'max';
end
switch lower(mode)
    case 'max'
        c = max(u);
    case 'none'
        c = 0;
    otherwise
        error('Unknown gauge mode "%s".', mode);
end
un = u - c;
scale = exp(min(max(c/eps, -700), 700));
Wn = W * scale;
end
