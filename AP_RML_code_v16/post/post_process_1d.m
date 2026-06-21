% POST_PROCESS_1D
% -------------------------------------------------------------------------
% 一维后处理统一入口。你可以在这里填一个或多个数据目录，也可以指定文件名。
% 画图逻辑由 plot_case_1d / plot_history_1d / plot_snapshot_1d 完成。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口：数据位置 =========================

% 写一个或多个数据目录。目录一般来自：
%   data/算例名/wkb_eps..._dt..._Nx..._Ntheta..._t...
%   data/算例名/direct_eps..._dt..._Nx..._Ntheta..._t...
dataDirs = { ...
    % fullfile(root, 'data', 'test3_wkb_advantage_refinement', 'wkb_eps0p01_dt0p001_Nx40_Ntheta16_t2'), ...
    % fullfile(root, 'data', 'test3_wkb_advantage_refinement', 'direct_eps0p01_dt0p001_Nx40_Ntheta32_t2') ...
    };

% 可选：指定每个目录下要画的文件名。
% 为空时自动使用 final 文件；若 opt.snapshotTime 有值，则自动找最接近该时刻的 snapshot。
fileNames = {};
% fileNames = {'snapshot_wkb_step0001000_t1.mat', 'snapshot_direct_step0001000_t1.mat'};

%% ========================= 用户接口：画什么、存到哪里 =========================

opt.jobs = {'history', 'snapshot'};  % 'history'、'snapshot'
opt.snapshotTime = [];               % [] 表示 final；例如 1 表示找 t≈1 的 snapshot
opt.outputDir = fullfile(root, 'data', 'figures', 'manual_post');
opt.figurePrefix = 'manual';
opt.savePng = true;
opt.saveFig = false;
opt.closeFigure = false;

%% ========================= 运行后处理 =========================

if isempty(dataDirs)
    fprintf(['请先在 post_process_1d.m 的 dataDirs 中填入一个或多个数据目录。\n' ...
             '数据目录格式通常为 data/算例名/wkb_eps... 或 data/算例名/direct_eps...。\n']);
else
    figInfo = plot_case_1d(dataDirs, fileNames, opt); %#ok<NASGU>
end
