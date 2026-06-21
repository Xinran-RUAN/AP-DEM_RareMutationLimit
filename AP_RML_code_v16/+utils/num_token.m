function s = num_token(x)
%NUM_TOKEN 将数值转换成不会丢精度的文件名片段。
%
%  例子：
%      0.01      -> 0p01
%      0.001     -> 0p001
%      2.5e-4    -> 2p5em4
%      1e-4      -> 1em4
%      2.5       -> 2p5
%      10000     -> 1ep4
%
%  旧版曾对小数使用 %.0e，这会把 2.5e-4 写成 3e-4。
%  这里统一保留足够有效数字，避免 eps、dt 等参数在文件名中被四舍五入。

if isempty(x) || ~isnumeric(x) || ~isscalar(x) || ~isfinite(x)
    s = 'nan';
    return;
end

x = double(x);

if x == 0
    s = '0';
    return;
end

ax = abs(x);

if ax < 1e-4 || ax >= 1e4
    raw = sprintf('%.12e', x);
    raw = regexprep(raw, '(\.\d*?)0+(e[+-]?\d+)$', '$1$2');
    raw = regexprep(raw, '\.(e[+-]?\d+)$', '$1');
else
    raw = sprintf('%.12g', x);
end

raw = lower(raw);
raw = strrep(raw, '+', '');
raw = regexprep(raw, 'e-0*', 'em');
raw = regexprep(raw, 'e0*', 'ep');
raw = strrep(raw, '.', 'p');
raw = strrep(raw, '-', 'm');
raw = regexprep(raw, '[^A-Za-z0-9]+', '_');
raw = regexprep(raw, '^_+|_+$', '');

if isempty(raw)
    raw = 'nan';
end

s = raw;
end
