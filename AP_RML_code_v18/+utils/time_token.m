function token = time_token(t)
%TIME_TOKEN 将时间转换成稳定、可读且不丢精度的文件名片段。
%
%  与 utils.num_token 使用同一套规则，例如
%      2.5   -> 2p5
%      0.5   -> 0p5
%      1e-4  -> 1em4
%  避免文件名中出现小数点，也避免 2.5e-4 被写成 3e-4。

token = utils.num_token(t);
end
