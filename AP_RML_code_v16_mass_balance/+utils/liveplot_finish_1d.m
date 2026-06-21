function liveplot_finish_1d(lp)
%LIVEPLOT_FINISH_1D 在线图像结束时的清理函数。
% 当前版本不主动关闭窗口，便于用户在计算结束后继续查看最终图像。
if nargin < 1 || isempty(lp)
    return;
end
if isfield(lp, 'enabled') && lp.enabled && isfield(lp, 'fig') && isgraphics(lp.fig)
    drawnow;
end
end
