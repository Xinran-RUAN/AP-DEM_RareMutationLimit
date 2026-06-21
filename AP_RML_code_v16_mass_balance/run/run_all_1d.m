% RUN_ALL_1D
% -------------------------------------------------------------------------
% 批量调用若干一维算例脚本。这里只列脚本名，不写数值参数。
% 每个脚本自己的“用户接口”仍在对应文件开头。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

scriptList = { ...
    'run_smoke_test_1d', ...
    'run_main_1d' ...
    % 'run_test1_assumption_diagnostics', ...
    % 'run_test2_long_time_concentration', ...
    % 'run_test3_wkb_advantage_refinement', ...
    % 'run_test4_coarse_theta_AP', ...
    % 'run_test5_accuracy' ...
    };

%% ========================= 运行 =========================

for i = 1:numel(scriptList)
    fprintf('\n========== [%d/%d] %s ==========%s', i, numel(scriptList), scriptList{i}, newline);
    run(scriptList{i});
end
