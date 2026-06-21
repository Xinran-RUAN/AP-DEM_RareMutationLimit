function str = format_seconds(seconds)
%FORMAT_SECONDS 把秒数格式化成便于阅读的时间字符串。
%   主要用于进度输出中的 elapsed / ETA 显示。
if nargin < 1 || isempty(seconds) || ~isfinite(seconds) || seconds < 0
    str = '--';
    return;
end
if seconds < 60
    str = sprintf('%.1fs', seconds);
elseif seconds < 3600
    m = floor(seconds/60);
    s = seconds - 60*m;
    str = sprintf('%dm%04.1fs', m, s);
else
    h = floor(seconds/3600);
    rems = seconds - 3600*h;
    m = floor(rems/60);
    s = rems - 60*m;
    str = sprintf('%dh%02dm%04.1fs', h, m, s);
end
end
