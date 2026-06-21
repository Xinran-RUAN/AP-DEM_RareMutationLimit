% RUN_TEST5_ACCURACY
% -------------------------------------------------------------------------
% 算例 5：空间网格加密精度检查。唯一批量变量是 Nx。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口 =========================

cfg.caseName = 'test5_accuracy';
cfg.solvers = {'wkb','direct'};

cfg.profile = 'baseline';
cfg.eps = 1e-2;
cfg.Nx = [20, 40, 80];       % 唯一批量变量
cfg.Ntheta_wkb = 64;
cfg.Ntheta_direct = 64;
cfg.dt = 1e-3;
cfg.T = 1;
cfg.tol = 1e-9;

cfg.u0_label = 'u_peak_070';
cfg.W0_label = 'W_flat';
cfg.u0_fun = @(theta) -0.05*(1 - cos(2*pi*(theta - 0.70)));
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

cfg.timeIntegrator = 'frozen';
cfg.phaseHamiltonian = 'weno5';
cfg.amplitudeVariant = 'split';
cfg.rhoReconstruction = 'laplace-hybrid';
cfg.reactionDiscretization = 'implicit';

cfg.snapshotTimes = [0, cfg.T];
cfg.historyEveryTime = 0.25;
cfg.clearOutput = true;
cfg.verbose = true;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
