function startup()
% Add project root to path
this = fileparts(mfilename('fullpath'));   % .../project_root/runs
root = fileparts(this);                    % .../project_root
addpath(root);                             % MUST add root (package folder lives here)
% 可选：如果你还有非 package 的普通函数目录，也可以加：
% addpath(genpath(root));
end