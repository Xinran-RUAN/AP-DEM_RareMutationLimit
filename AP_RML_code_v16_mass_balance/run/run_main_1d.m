% RUN_MAIN_1D
% -------------------------------------------------------------------------
% 普通一维主程序。所有人为可调参数都集中在本文件开头。
% 真正运行由 run_case_1d 完成；本脚本只负责描述一个算例。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口：算例名称和求解器 =========================

cfg.caseName = 'main_1d';
cfg.solvers = {'wkb'};        % 可选：{'wkb'}、{'direct'}、{'wkb','direct'}

% 若某个参数需要批量跑，就把它写成列表。一个主程序最多只批量跑一个变量。
% 例：cfg.eps = [1e-2, 1e-3]; 或 cfg.Ntheta_direct = [32, 64, 128];

%% ========================= 用户接口：基本数值参数 =========================

cfg.profile = 'baseline';
cfg.eps = 1e-2;
cfg.Nx  = 50;
cfg.Ntheta_wkb = 32;
cfg.Ntheta_direct = 64;
cfg.dt  = 1e-3;
cfg.T   = 2;
cfg.tol = 1e-8;

%% ========================= 用户接口：模型系数 =========================

cfg.K_label = 'K1';
cfg.D_label = 'D1';

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

% 平坦相位初值示例：
% cfg.u0_label = 'u_flat';
% cfg.u0_fun = @(theta) 0.*theta;

%% ========================= 用户接口：算法选择 =========================

cfg.timeIntegrator = 'frozen';            % 'frozen', 'projective', 'implicit-coupled'
cfg.phaseHamiltonian = 'weno5';           % 'weno5', 'godunov', 'lf'
cfg.amplitudeVariant = 'split';           % 'split', 'wkb-if-split', 'density-compatible'
cfg.rhoReconstruction = 'laplace-hybrid'; % 'direct-log', 'laplace', 'laplace-hybrid'
cfg.reactionDiscretization = 'implicit';  % 'implicit', 'logistic-exact'
cfg.residualMode = 'legacy';              % 'legacy', 'wkb-steady'

% 质量平衡修正：每步用 eps*dM/dt = int rho*(K-rho) dx 修正整体质量。
cfg.massCorrectionDuringSolve = true;       % {true,false}
cfg.massCorrectionFormula = 'trapezoid';    % {'trapezoid','backward-euler'}

cfg.adaptiveTimeStep = false;
cfg.stopByResidual = false;

%% ========================= 用户接口：保存和输出 =========================

cfg.snapshotTimes = unique([0, cfg.T]);
cfg.snapshotEveryTime = Inf;
cfg.saveLatestCheckpoint = true;
cfg.checkpointEveryTime = Inf;
cfg.saveHistoryCheckpoint = true;
cfg.historyEveryTime = 0.5;
cfg.saveFinalResult = true;

cfg.clearOutput = true;
cfg.verbose = true;
cfg.progressEveryPercent = 5;
cfg.progressEverySeconds = 20;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
