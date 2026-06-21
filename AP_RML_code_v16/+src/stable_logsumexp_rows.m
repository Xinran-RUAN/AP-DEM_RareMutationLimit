function [sumVal, logSumVal, info] = stable_logsumexp_rows(logA, weights, logCap)
%STABLE_LOGSUMEXP_ROWS  对每一行做带权 log-sum-exp。
%
%  sumVal_j = sum_k weights_k * exp(logA_{j,k})
%  该函数避免中间 exp 溢出。若真实 logSum 超过 logCap，会把输出截断到
%  exp(logCap)，并在 info.numCappedHigh 中记录。截断只用于保护诊断/重构，
%  真正发生截断时应查看 history 中的监测量，而不是把结果当作可靠收敛。

if nargin < 2 || isempty(weights)
    weights = ones(1,size(logA,2));
end
if nargin < 3 || isempty(logCap)
    logCap = 690;
end
weights = weights(:).';
logW = log(max(weights, realmin));
logB = bsxfun(@plus, real(logA), logW);
rowMax = max(logB, [], 2);
rowMaxSafe = rowMax;
rowMaxSafe(~isfinite(rowMaxSafe)) = -Inf;
scaled = exp(bsxfun(@minus, logB, rowMaxSafe));
scaled(~isfinite(scaled)) = 0;
s = sum(scaled, 2);
logSumVal = rowMaxSafe + log(max(s, realmin));

info = struct();
info.logMaxInput = max(logA(:));
info.logMinInput = min(logA(:));
info.logSumMax = max(logSumVal(:));
info.logSumMin = min(logSumVal(:));
info.numCappedHigh = nnz(logSumVal > logCap);
info.numNonfiniteLogInput = nnz(~isfinite(logA(:)));

logOut = min(logSumVal, logCap);
logOut(~isfinite(logOut)) = -Inf;
sumVal = exp(logOut);
end
