function B = interp_periodic_theta_1d(thetaSrc, A, thetaDst, method)
%INTERP_PERIODIC_THETA_1D 周期 theta 插值。
% A 可以是 1 x K 或 nx x K，返回尺寸为对应的 1 x M 或 nx x M。
%
% 说明：
%   MATLAB 的 interp1('spline'/'pchip') 在某些退化输入下会调用
%   griddedInterpolant 并报错：
%       插值要求每个网格维度至少有两个采样点。
%   双网格长时间推进中，一旦某个中间量退化为单点/空点，这个错误会
%   中断主循环。这里做了防护：若有效 theta 点少于方法要求，自动退化
%   为 linear/constant；若 spline/pchip/makima 失败，也自动 fallback。

if nargin < 4 || isempty(method)
    method = 'spline';
end
method = lower(char(method));

if nargin < 3 || isempty(thetaDst)
    % 保持行数，返回空 theta 维。
    A = real(A);
    if isvector(A)
        B = zeros(1,0);
    else
        B = zeros(size(A,1),0);
    end
    return;
end

if isempty(thetaSrc)
    error('interp_periodic_theta_1d:EmptySourceGrid', 'thetaSrc 为空，无法插值。');
end

thetaSrc = mod(thetaSrc(:), 1);
thetaDst = mod(thetaDst(:), 1);
A = real(A);
wasVector = isvector(A);
if wasVector
    Arow = A(:).';
else
    Arow = A;
end

% 允许 A 为 K x nx 的转置形式。
if size(Arow,2) ~= numel(thetaSrc)
    if size(Arow,1) == numel(thetaSrc)
        Arow = Arow.';
    else
        error('interp_periodic_theta_1d:SizeMismatch', ...
            'A 的 theta 维度与 thetaSrc 不匹配：size(A)=[%s], numel(thetaSrc)=%d。', ...
            num2str(size(A)), numel(thetaSrc));
    end
end

% 去掉非有限网格点。
goodTheta = isfinite(thetaSrc);
thetaSrc = thetaSrc(goodTheta);
Arow = Arow(:, goodTheta);
if isempty(thetaSrc)
    error('interp_periodic_theta_1d:NoFiniteSourceGrid', 'thetaSrc 没有有限采样点。');
end

[thetaSrc, idx] = sort(thetaSrc, 'ascend');
Arow = Arow(:, idx);
[thetaSrc, uniq] = unique(thetaSrc, 'stable');
Arow = Arow(:, uniq);

K = numel(thetaSrc);
M = numel(thetaDst);
if K <= 1
    % 单个采样点只能做常数延拓。
    B = repmat(Arow(:,1), 1, M);
else
    thetaExt = [thetaSrc(end)-1; thetaSrc; thetaSrc(1)+1];
    AExt = [Arow(:,end), Arow, Arow(:,1)];

    % 对退化小样本，避免 spline/makima 的阶数要求。
    methodTry = method;
    if K < 3 && any(strcmp(methodTry, {'spline','makima','pchip','cubic'}))
        methodTry = 'linear';
    elseif K < 4 && any(strcmp(methodTry, {'spline','makima','cubic'}))
        methodTry = 'pchip';
    end

    B = local_interp_with_fallback(thetaExt, AExt, thetaDst, methodTry);
end

B = real(B);
if wasVector
    B = B(:).';
end
end

function B = local_interp_with_fallback(thetaExt, AExt, thetaDst, methodTry)
% 返回 rows x numel(thetaDst)。
methods = {methodTry, 'pchip', 'linear', 'nearest'};
% 去重，保持顺序。
methods = methods(~cellfun(@isempty, methods));
[~, ia] = unique(methods, 'stable');
methods = methods(sort(ia));

lastME = [];
for im = 1:numel(methods)
    m = methods{im};
    try
        Bt = interp1(thetaExt, AExt.', thetaDst, m, 'extrap');
        B = Bt.';
        return;
    catch ME
        lastME = ME; %#ok<NASGU>
    end
end

% 最后防线：不用 interp1，直接常数延拓，避免长时间计算中断。
% 取原始周期块的第一个内点作为常数值。
try
    A0 = AExt(:,2);
catch
    A0 = AExt(:,1);
end
B = repmat(A0, 1, numel(thetaDst));
end
