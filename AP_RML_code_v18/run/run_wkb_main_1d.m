% RUN_WKB_MAIN_1D
% -------------------------------------------------------------------------
% 一维 WKB 主程序。
%
% 重要约定：
%   本代码不再使用 cfg.profile 或任何参数预设文件。
%   所有模型参数、初值、网格、时间步、算法选项、保存选项都在本文件中显式给出。
%
% WKB 双 theta 网格时：
%   cfg.Ntheta_wkb 是 W 的粗 theta 网格点数；
%   cfg.Ntheta_u   是 u 的细 theta 网格点数。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口：算例和输出根目录 =========================

cfg.caseName = 'main_1d';
cfg.solvers = {'wkb'};                    % 固定：本主程序只跑 WKB
cfg.dataRoot = fullfile(root, 'data');

%% ========================= 用户接口：基本数值参数 =========================

cfg.eps = 1e-5;
cfg.Nx  = 50;

% WKB 双 theta 网格：W 粗网格 + u 细网格。
cfg.useDualThetaGrid = true;
cfg.Ntheta_wkb = 16;       % W 的粗 theta 网格点数
cfg.Ntheta_u   = 64;       % u 的细 theta 网格点数
cfg.Ntheta_W   = [];       % [] 表示使用 cfg.Ntheta_wkb 作为 W 粗网格

cfg.dt  = 1e-2;
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

cfg.u0_label = 'u_flat';
cfg.W0_label = 'W_flat';

cfg.u0_fun = @(theta) 0.*theta;
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

% 非平坦相位初值示例：
% cfg.u0_label = 'u_peak_070';
% cfg.u0_fun = @(theta) -0.05*(1 - cos(2*pi*(theta - 0.70)));

% 初始 gauge 和初始质量处理。
cfg.initialGaugeMode = 'max';             % {'max','none'}
cfg.prepareInitialMass = false;           % {true,false}
cfg.initialRhoProfile = 'none';           % 仅 prepareInitialMass=true 时作为标签保留
cfg.initialRhoScale = 1.0;

%% ========================= 用户接口：WKB 时间推进和相位格式 =========================

cfg.timeIntegrator = 'frozen';             % 当前双网格主路径使用 frozen-H 时间推进
cfg.phaseHamiltonian = 'weno5';            % {'weno5','godunov','lf'}
cfg.lfAlpha = 1.0;                         % phaseHamiltonian='lf' 时使用；其它模式保留为显式参数

%% ========================= 用户接口：W 振幅更新 =========================

% full-eigen-relax：
%   保留完整 W 更新格式，包括 x 扩散、theta 扩散、theta 输运和反应项。
%   只是 W 粗网格上的 H 直接在 W 网格上求，避免从 u 网格插值 H 破坏
%   小 eps 时的 x 方向特征松弛结构。
cfg.amplitudeVariant = 'full-eigen-relax'; % {'split','full-eigen-relax'}
cfg.transportImplicit = true;              % 完整 W 格式中的 theta 输运项隐式
cfg.clipNegativeW = true;                  % 线性求解后把极小负值裁剪到 WFloor
cfg.WFloor = 0;                            % W 的非负下界
cfg.monitorMatrixCondition = false;        % 是否估计 W 线性系统条件数，较慢

%% ========================= 用户接口：rho 重构 =========================

cfg.rhoReconstruction = 'laplace-hybrid';  % {'direct-log','laplace','laplace-hybrid'}

%% ========================= 用户接口：反应项和残差 =========================

cfg.reactionDiscretization = 'implicit';   % {'implicit','logistic-exact'}
cfg.residualMode = 'legacy';               % {'legacy','wkb-steady'}

%% ========================= 用户接口：质量修正 =========================

cfg.massCorrectionDuringSolve = true;      % {true,false}
cfg.massCorrectionFormula = 'trapezoid';   % {'trapezoid','backward-euler'}
cfg.massCorrectionMinScale = 0;
cfg.massCorrectionMaxScale = Inf;
cfg.massCorrectionVerbose = false;

%% ========================= 用户接口：Hamiltonian / 特征值求解 =========================

% u 细网格上的 H 求解方式。
cfg.HSolverMode = 'direct-unique-D';       % {'direct-all-theta','direct-unique-D','D-pchip','D-pchip-local-correct'}
cfg.HDtol = 1e-12;
cfg.HInterpND = 24;
cfg.HLocalRadius = 4;
cfg.eigenSolver = 'auto';                  % {'auto','eig-full','eigs'}
cfg.eigenDenseThreshold = 180;
cfg.eigenTol = 1e-10;
cfg.eigenMaxit = 300;

% W 粗网格上的 H：full-eigen-relax 分支会直接在 W 网格求。
cfg.dualHSolverModeW = 'direct-unique-D';  % {'direct-all-theta','direct-unique-D'}
cfg.hamiltonianGauge = 'none';             % full-eigen-relax 建议先用 none
cfg.fullEigenRelaxHMode = 'raw';           % {'raw','consistent-gauge'}

%% ========================= 用户接口：双 theta 网格插值 =========================

cfg.dualWtoUInterp = 'spline';             % {'spline','pchip','linear','makima'}
cfg.dualUtoWInterp = 'spline';             % {'spline','pchip','linear','makima'}
cfg.dualHtoWInterp = 'spline';             % {'spline','pchip','linear','makima'}

%% ========================= 用户接口：时间步控制和终止条件 =========================

cfg.adaptiveTimeStep = false;
cfg.adaptiveStrategy = 'none';
cfg.adaptiveSafety = 0.8;
cfg.phaseCflSafety = 0.8;
cfg.dtMax = cfg.dt;
cfg.dtMin = 0;
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
