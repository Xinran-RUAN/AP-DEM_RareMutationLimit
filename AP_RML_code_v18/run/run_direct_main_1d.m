% RUN_DIRECT_MAIN_1D
% -------------------------------------------------------------------------
% 一维 direct density 主程序。
%
% 重要约定：
%   本代码不再使用 cfg.profile 或任何参数预设文件。
%   所有模型参数、初值、网格、时间步、算法选项、保存选项都在本文件中显式给出。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口：算例和输出根目录 =========================

cfg.caseName = 'main_1d';
cfg.solvers = {'direct'};                 % 固定：本主程序只跑 direct
cfg.dataRoot = fullfile(root, 'data');

%% ========================= 用户接口：基本数值参数 =========================

cfg.eps = 1e-5;
cfg.Nx  = 50;
cfg.Ntheta_direct = 64;
cfg.dt  = 1e-3;
cfg.T   = 10;
cfg.tol = 1e-8;
cfg.maxSteps = Inf;

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

% direct 初值通过同一套 WKB 型初值构造：
%     n0(x,theta) = W0(x,theta) * exp(u0(theta)/eps)。
% 这样 direct 与 WKB 可以使用完全一致的初值。
cfg.u0_label = 'u_flat';
cfg.W0_label = 'W_flat';

cfg.u0_fun = @(theta) 0.*theta;
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

% 初始 gauge 和初始质量处理。
cfg.initialGaugeMode = 'max';             % {'max','none'}
cfg.prepareInitialMass = false;           % {true,false}
cfg.initialRhoProfile = 'none';
cfg.initialRhoScale = 1.0;

%% ========================= 用户接口：direct 算法选择 =========================

cfg.reactionDiscretization = 'implicit';  % {'implicit','patankar','logistic-exact'}
cfg.residualMode = 'legacy';              % {'legacy'}

%% ========================= 用户接口：质量修正 =========================

cfg.massCorrectionDuringSolve = true;      % {true,false}
cfg.massCorrectionFormula = 'trapezoid';   % {'trapezoid','backward-euler'}
cfg.massCorrectionMinScale = 0;
cfg.massCorrectionMaxScale = Inf;
cfg.massCorrectionVerbose = false;

%% ========================= 用户接口：时间步控制和终止条件 =========================

cfg.adaptiveTimeStep = false;
cfg.stopByResidual = false;
cfg.progressPrintWhenTolReached = false;

%% ========================= 用户接口：保存和输出 =========================

cfg.storeSnapshots = true;
cfg.saveIntermediate = true;
cfg.snapshotTimes = unique([0, cfg.T]);
cfg.snapshotEveryTime = Inf;
cfg.snapshotEverySteps = [];
cfg.saveInitialSnapshot = true;
cfg.saveInitialAcceptedSnapshot = false;
cfg.saveLatestCheckpoint = true;
cfg.checkpointEveryTime = Inf;
cfg.checkpointEverySteps = [];
cfg.saveHistoryCheckpoint = true;
cfg.historyEveryTime = 0.5;
cfg.historyEverySteps = [];
cfg.saveFinalResult = true;
cfg.saveMatV73 = false;

cfg.clearOutput = true;
cfg.continueOnError = false;
cfg.keepResults = false;
cfg.dryRun = false;

cfg.verbose = true;
cfg.progressEveryPercent = 5;
cfg.progressEverySeconds = 20;
cfg.progressEverySteps = [];
cfg.printIntermediateSaveMessage = false;
cfg.livePlot = false;
cfg.livePlotEveryTime = 0.5;
cfg.livePlotEverySteps = [];
cfg.livePlotSave = false;
cfg.livePlotPause = 0.0;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>
