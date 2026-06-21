% RUN_TEMPLATE_1D
% -------------------------------------------------------------------------
% 新算例模板。复制本文件后，只改“用户接口”部分即可。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口：算例名称和求解器 =========================

cfg.caseName = 'template_1d';
cfg.solvers = {'wkb','direct'};

%% ========================= 用户接口：基本数值参数 =========================

cfg.profile = 'baseline';
cfg.eps = 1e-2;
cfg.Nx  = 40;
cfg.Ntheta_wkb = 32;
cfg.Ntheta_direct = 64;       % 可写成列表，例如 [64, 128]
cfg.dt  = 1e-3;
cfg.T   = 2;
cfg.tol = 1e-8;

%% ========================= 用户接口：模型系数 =========================

cfg.K_label = 'K_default';
cfg.D_label = 'D_default';

cfg.K0 = 1.0;
cfg.K1 = 0.5;
cfg.Dmin = 0.2;
cfg.D1 = 0.4;
cfg.theta_m = 0.35;

K0 = cfg.K0;
K1 = cfg.K1;
Dmin = cfg.Dmin;
D1 = cfg.D1;
theta_m = cfg.theta_m;

cfg.K_fun = @(x) K0 - K1*cos(2*pi*x);
cfg.D_fun = @(theta) Dmin + D1*(1 - cos(2*pi*(theta - theta_m)));

%% ========================= 用户接口：初值 =========================

cfg.u0_label = 'u_peak_070';
cfg.W0_label = 'W_flat';

cfg.u0_fun = @(theta) -0.05*(1 - cos(2*pi*(theta - 0.70)));
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

%% ========================= 用户接口：算法选择 =========================

cfg.timeIntegrator = 'frozen';
cfg.phaseHamiltonian = 'weno5';
cfg.amplitudeVariant = 'split';
cfg.rhoReconstruction = 'laplace-hybrid';
cfg.reactionDiscretization = 'implicit';
cfg.residualMode = 'legacy';

cfg.adaptiveTimeStep = false;
cfg.stopByResidual = false;

%% ========================= 用户接口：保存和输出 =========================

cfg.snapshotTimes = [0, 0.5, 1, cfg.T];
cfg.snapshotTimes = unique(cfg.snapshotTimes(cfg.snapshotTimes <= cfg.T + 1e-12));
cfg.saveFinalResult = true;
cfg.clearOutput = true;
cfg.verbose = true;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
