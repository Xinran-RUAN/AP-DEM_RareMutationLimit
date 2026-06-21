function [P, nFine, WFine, UFine, thetaFine] = reconstruct_trait_marginal_from_wkb_1d(W, u, op, eps, thetaFine, method)
%RECONSTRUCT_TRAIT_MARGINAL_FROM_WKB_1D 在细 theta 网格上由 W,u 正确重构 P(theta)。
%
%  这是为后处理画图专门保留的安全接口。正确顺序必须是：
%      1. 分别把 W(x,theta) 和 u(theta) 插值到目标 thetaFine；
%      2. 在 thetaFine 上重构 nFine = WFine .* exp(UFine/eps)；
%      3. 再对 x 积分得到 P(thetaFine)。
%
%  不要先在粗网格上算 P_k = int_x W exp(u/eps) dx 再插值 P_k。
%  对小 epsilon，P(theta) 可能非常尖，直接插值 P 会把“指数相位”与
%  “振幅”混在一起，既会错峰，也会错质量；旧版 post 脚本中的做法是
%  本函数采用的做法。
%
%  输入：
%    W          : [Nx+1, K] 振幅；
%    u          : [1,K] 或 [Nx+1,K] 相位；
%    op         : build_operators_1d 返回的结构，含 x/theta/dx；
%    eps        : 小突变参数；
%    thetaFine  : 目标 theta 网格；为空时默认取 max(8*K,512) 个点；
%    method     : 'spline','pchip','linear' 等 interp1 方法，默认 'spline'。
%
%  输出：
%    P          : 目标 theta 网格上的 trait marginal；
%    nFine      : 细 theta 网格上的重构密度；
%    WFine,UFine: 插值后的 W,u；
%    thetaFine  : 实际使用的目标网格。

if nargin < 5 || isempty(thetaFine)
    Kfine = max(8*op.Ntheta, 512);
    thetaFine = linspace(0, 1, Kfine+1).';
    thetaFine(end) = [];
else
    thetaFine = thetaFine(:);
end
if nargin < 6 || isempty(method)
    method = 'spline';
end

theta = op.theta(:);
xn = size(W,1);
kn = numel(thetaFine);
WFine = zeros(xn, kn);
for j = 1:xn
    WFine(j,:) = local_periodic_interp(theta, W(j,:).', thetaFine, 1.0, method).';
end
WFine(~isfinite(WFine)) = 0;
WFine = max(real(WFine), 0);

if isvector(u)
    uFineRow = local_periodic_interp(theta, u(:), thetaFine, 1.0, method).';
    UFine = repmat(uFineRow, xn, 1);
else
    UFine = zeros(xn, kn);
    for j = 1:xn
        UFine(j,:) = local_periodic_interp(theta, u(j,:).', thetaFine, 1.0, method).';
    end
end
UFine = real(UFine);

arg = UFine / eps;
arg = min(max(arg, -745), 700);  % 双精度 exp 安全范围
nFine = WFine .* exp(arg);
nFine(~isfinite(nFine)) = 0;

wx = utils.trapz_weights(op.nx, op.dx);
P = wx.' * nFine;
P = real(P(:));
end

function fq = local_periodic_interp(theta, f, thetaq, period, method)
% 周期插值：在末端补一个周期点，避免 theta=0/1 边界断裂。
theta = theta(:);
f = f(:);
thetaq = mod(thetaq(:), period);
[theta, idx] = sort(mod(theta, period));
f = f(idx);

% 去掉可能重复的周期端点。
if numel(theta) > 1 && abs(theta(end)-period) < 1e-12
    theta = theta(1:end-1);
    f = f(1:end-1);
end
[theta, uniq] = unique(theta, 'stable');
f = f(uniq);

thetaExt = [theta; theta(1)+period];
fExt = [f; f(1)];
try
    fq = interp1(thetaExt, fExt, thetaq, method);
catch
    fq = interp1(thetaExt, fExt, thetaq, 'linear');
end
end
