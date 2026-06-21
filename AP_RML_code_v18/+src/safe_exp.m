function y = safe_exp(x, clipValue)
%SAFE_EXP Exponential with optional clipping to avoid Inf in diagnostic ratios.
%   y = SAFE_EXP(x) evaluates exp(x) after clipping x to [-700,700].
%   y = SAFE_EXP(x, clipValue) uses a user-specified clipping value.
%
%   The clipping is only a numerical safeguard against floating-point overflow.
%   In well-resolved WKB runs the exponents entering density-compatible ratios
%   should stay far from this bound.
if nargin < 2 || isempty(clipValue)
    clipValue = 700;
end
if ~isfinite(clipValue) || clipValue <= 0
    y = exp(x);
    return;
end
y = exp(min(max(x, -clipValue), clipValue));
end
