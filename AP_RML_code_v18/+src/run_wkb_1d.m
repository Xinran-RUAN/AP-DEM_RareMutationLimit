function result = run_wkb_1d(par)
%RUN_WKB_1D 运行一维 WKB 求解器。
%
% 本整理版只保留当前文章使用的双 theta 网格 WKB 主路径：
%   u(theta,t) 在细 theta 网格；
%   W(x,theta,t) 在粗 theta 网格；
%   W 更新可使用 full-eigen-relax 完整格式。
%
% 旧版 single-grid/time-AP/projective 补丁路径已从主入口移除，避免参数隐藏在
% 预设或自动补丁中。

if nargin < 1 || ~isstruct(par)
    error('run_wkb_1d:InvalidInput', 'par 必须是 struct。');
end
if ~isfield(par, 'useDualThetaGrid') || isempty(par.useDualThetaGrid) || ~logical(par.useDualThetaGrid)
    error('run_wkb_1d:DualThetaRequired', ...
        ['本整理版 WKB 主程序只保留双 theta 网格。请在 run_wkb_main_1d.m 中设置 ', ...
         'cfg.useDualThetaGrid=true, cfg.Ntheta_wkb 和 cfg.Ntheta_u。']);
end

result = src.run_wkb_dual_1d(par);
end
