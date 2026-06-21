function root = startup_AP_RML(~)
%STARTUP_AP_RML Add the folders needed by the 1D AP-RML code.
% 本整理版只保留 1D WKB、1D direct 和统一后处理；旧 2D/legacy 路径已删除。
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'run'));
addpath(fullfile(root, 'post'));
end
