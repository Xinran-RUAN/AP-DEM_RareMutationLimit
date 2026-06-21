% RUN_SMOKE_TEST_1D
% -------------------------------------------------------------------------
% 很小的冒烟测试：快速检查 wkb/direct 是否能跑通和保存文件。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

cfg.caseName = 'smoke_test_1d';
cfg.solvers = {'wkb','direct'};

cfg.profile = 'baseline';
cfg.eps = 1e-2;
cfg.Nx = 10;
cfg.Ntheta_wkb = 12;
cfg.Ntheta_direct = 12;
cfg.dt = 2e-3;
cfg.T = 0.02;
cfg.tol = 1e-8;

cfg.u0_label = 'u_flat';
cfg.W0_label = 'W_flat';
cfg.u0_fun = @(theta) 0.*theta;
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

cfg.timeIntegrator = 'frozen';
cfg.phaseHamiltonian = 'weno5';
cfg.amplitudeVariant = 'split';
cfg.rhoReconstruction = 'laplace-hybrid';
cfg.reactionDiscretization = 'implicit';

cfg.snapshotTimes = [0, cfg.T];
cfg.historyEveryTime = cfg.T;
cfg.checkpointEveryTime = Inf;
cfg.clearOutput = true;
cfg.verbose = true;
cfg.progressEveryPercent = 50;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
