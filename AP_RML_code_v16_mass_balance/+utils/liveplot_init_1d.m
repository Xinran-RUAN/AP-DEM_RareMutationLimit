function lp = liveplot_init_1d(par, mode)
%LIVEPLOT_INIT_1D 初始化一维计算过程中的动态图窗口。
%   mode = 'wkb' 或 'direct'。
%
% 说明：
%   1. 该函数仅在 par.livePlot=true 时被调用；
%   2. 如果当前 MATLAB 环境不支持图形窗口，则自动关闭 live plot；
%   3. 为避免拖慢计算，建议配合较大的 livePlotEveryTime 或
%      livePlotEverySteps 使用，而不是每一步都刷新。

if nargin < 2 || isempty(mode)
    mode = 'wkb';
end

lp = struct();
lp.enabled = false;
lp.mode = lower(mode);
lp.fig = [];
lp.frameDir = '';
lp.lastStep = 0;
lp.lastTime = 0;
lp.counter = 0;

if ~isfield(par, 'livePlot') || ~par.livePlot
    return;
end

try
    % 为在线图像单独建立保存目录。
    if isfield(par, 'livePlotSave') && par.livePlotSave
        lp.frameDir = utils.ensure_dir(fullfile(par.outdir, 'live_frames'));
    end

    switch lp.mode
        case 'wkb'
            figName = 'AP-RML 1D WKB live monitor';
        otherwise
            figName = 'AP-RML 1D direct live monitor';
    end

    lp.fig = figure('Name', figName, 'NumberTitle', 'off', ...
        'Color', 'w', 'Visible', 'on');
    clf(lp.fig);
    lp.enabled = true;
catch ME
    warning('liveplot:initFailed', ...
        '无法初始化在线图像窗口，已自动关闭 livePlot。原因: %s', ME.message);
    lp.enabled = false;
end
end
