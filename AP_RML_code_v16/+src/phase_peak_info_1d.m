function info = phase_peak_info_1d(u, op)
%PHASE_PEAK_INFO_1D  用三点二次插值估计相位 u 的连续峰值。
%
%  WKB 方法关心的是相位最大点。小 epsilon 时，密度峰宽可能远小于
%  trait 网格间距，单纯取 max(u_k) 会低估真正的连续峰值，并使
%  exp(u/eps) 的重构对网格点位置非常敏感。本函数在 nodal 最大点附近
%  使用周期三点二次插值，返回：
%     k0         : nodal 最大点编号；
%     thetaStar  : 二次插值得到的亚网格峰值位置；
%     zStar      : 相对于 theta(k0) 的局部坐标，范围限制在 [-dtheta,dtheta]；
%     uStar      : 二次插值得到的峰值；
%     dduStar    : 二阶导数估计；
%     ok         : 是否为凹峰。
%
%  该函数只用于 gauge、Laplace 重构和诊断，不改变主格式的 WKB 变量。

K = op.Ntheta;
dth = op.dtheta;
u = real(u(:).');

[~, k0] = max(u);
km = mod(k0-2, K) + 1;
kp = mod(k0,   K) + 1;

um = u(km); u0 = u(k0); up = u(kp);
denom = up - 2*u0 + um;
ddu = denom / dth^2;

ok = isfinite(ddu) && ddu < 0 && isfinite(um) && isfinite(u0) && isfinite(up);
if ok && abs(denom) > realmin
    zStar = -0.5*dth*(up - um)/denom;
    % 为了避免在非单峰或噪声情况下跑到相邻单元外，把顶点限制在局部 stencil 内。
    zStar = max(min(zStar, dth), -dth);
else
    zStar = 0;
    ok = false;
end

% Lagrange 插值权重，对应局部节点 {-dth,0,dth}。
zm = -dth; z0 = 0; zp = dth;
L1 = (zStar - z0)*(zStar - zp)/((zm-z0)*(zm-zp));
L2 = (zStar - zm)*(zStar - zp)/((z0-zm)*(z0-zp));
L3 = (zStar - zm)*(zStar - z0)/((zp-zm)*(zp-z0));
uStar = um*L1 + u0*L2 + up*L3;

info = struct();
info.k0 = k0;
info.km = km;
info.kp = kp;
info.thetaNode = op.theta(k0);
info.thetaStar = mod(op.theta(k0) + zStar, 1);
info.zStar = zStar;
info.uNodeMax = u0;
info.uStar = uStar;
info.dduStar = ddu;
info.ok = ok;
info.L = [L1, L2, L3];
end
