function par = default_params_1d(profile)
%DEFAULT_PARAMS_1D 一维 WKB/direct 程序的默认参数。
%
%  本合并版的设计原则：
%    1. 默认完全偏向“旧版/论文基准”行为，便于复现实验和旧脚本；
%    2. v13 中加入的小 epsilon、时间 AP、checkpoint 等补丁全部保留，
%       但默认作为可选项关闭或不自动改写用户设置；
%    3. 若要尝试 v13 实验补丁，建议调用
%           par = model.apply_patch_preset_1d(par,'v13-small-eps');
%       或在 run_main_1d.m 中把 experimentPreset 改为 'v13-small-eps'、
%       'v13-projective'、'v13-implicit'。
%
%  profile = 'baseline'       旧版基准模型，K(x)=K0-K1*cos(2*pi*x)；
%  profile = 'baseline-v13'   v13 曾使用的 K(x)=K0+K1*cos(2*pi*x)，仅作对比；
%  profile = 'heterogeneous'  非平凡非均匀算例。

if nargin < 1 || isempty(profile)
    profile = 'baseline';
end
par.profile = profile;

%=========================================================================
% 1. 网格和时间推进参数
%=========================================================================
par.Nx = 80;                  % x 方向子区间数；实际节点数为 Nx+1
par.Ntheta = 64;              % 周期 trait 网格点数
par.dt = 1e-3;                % 固定时间步；若启用自适应时间步，它作为最大时间步
par.T = 20;                   % 终止时间
par.maxSteps = Inf;           % 额外最大步数；Inf 表示只受 T 控制
par.tol = 1e-8;               % 残差停止阈值
par.saveEvery = Inf;          % 预留字段
par.outdir = utils.output_dir(profile);

%=========================================================================
% 2. 小突变参数和初值参数
%=========================================================================
par.eps = 1e-2;
par.theta0 = 0.70;            % peaked 初值时的相位峰值位置；旧版 flat 初值不使用
par.alpha0 = 0.05;            % peaked 初值的相位宽度参数

% 旧版当前副本的有效初值是 flat：u(theta)=0，W(x,theta)=1。
% 这不会预先指定 trait 峰，适合 long-time concentration 实验。
% 如需测试原始 peaked WKB 初值，可设为 'peaked'。
par.initialPhaseMode = 'flat';        % 'flat'/'homogeneous' 或 'peaked'
par.initialAmplitudeMode = 'flat-rho';% 'flat-rho': W=1；'mild-x': W=1+0.1*cos(2*pi*x)

% 初始 gauge 默认采用旧版 nodal max 归一化。v13 的 phase-peak gauge
% 对小 epsilon 更稳，但会略改初始 W/u 的分配，因此默认不启用。
par.initialGaugeMode = 'legacy-max';  % 'legacy-max', 'phase-peak', 'none'

% v13 的 well-prepared mass 补丁：只缩放 W，使初始 rho 为 O(1)。
% 默认关闭，避免自动改变旧版初值；需要小 epsilon 稳健实验时再打开。
par.prepareInitialMass = false;
par.initialRhoProfile = 'scaledK';    % prepareInitialMass=true 时：'scaledK','constant','current'
par.initialRhoScale = 0.5;

%=========================================================================
% 3. 时间积分器实验选项（v13/v11 补丁，默认关闭）
%=========================================================================
% 'frozen'            : 旧版路径，H=H(rho^n) 在一个时间步内冻结；
% 'implicit-coupled'  : 实验版全隐式耦合 Picard；
% 'projective'        : 实验版 micro-macro/projective。
par.timeIntegrator = 'frozen';

% 方案 1：全隐式耦合 Picard 参数。
par.implicitMaxIter = 12;
par.implicitTol = 1e-7;
par.implicitRelax = 0.7;
par.implicitMinIter = 2;
par.implicitRetryOnFail = true;
par.implicitFailureTol = 5e-2;

% 方案 3：projective / micro-macro 参数。
par.projectiveMicroSteps = 8;
par.projectiveMicroDtFactor = 0.25;
par.projectiveMinMicroDt = 1e-12;
par.projectiveMaxMicroDt = Inf;
par.projectiveMode = 'u-only';
par.projectiveDerivativeClipU = 50;
par.projectiveDerivativeClipLogW = 50;
par.projectiveLogWAbsClip = 650;
par.projectiveCorrectorSteps = 0;
par.projectiveUseMacroRhoRefresh = true;

%=========================================================================
% 4. WKB 数值格式选择
%=========================================================================
par.phaseHamiltonian = 'weno5';       % 旧版默认：'godunov', 'lf', 或 'weno5'
par.lfAlpha = 10;

% 旧版默认 split。density-compatible 现在明确映射到旧半离散矩阵版本；
% 若要用 v13 的 integrating factor WKB 格式，请显式设为 'wkb-if-split'。
par.amplitudeVariant = 'split';       % 'split','density-compatible','wkb-if-split','split-if'

% direct solver 的旧版反应离散等价于 linear implicit。WKB split 旧格式不使用
% reactionDiscretization；wkb-if-split 新格式推荐改为 'logistic-exact'。
par.reactionDiscretization = 'implicit';
par.transportImplicit = true;
par.rhoReconstruction = 'laplace-hybrid'; % 旧版默认：'direct-log','laplace','laplace-hybrid'
par.residualMode = 'legacy';   % 'legacy' 或 'wkb-steady'

% 诊断中 theta_wkb 的定义。旧版用节点最大值；v13 可用二次插值亚网格峰。
par.useSubgridPhasePeak = false;

%=========================================================================
% 5. 小 epsilon 稳健化补丁（全部可选）
%=========================================================================
% 主开关：false 时，run_wkb_1d 不会因为 eps 小而自动把 direct-log 改为
% laplace-hybrid，也不会自动把 reaction/gauge 改成 v13 推荐值。
par.smallEpsAutoPatches = false;
par.enableV13AutoPatches = false;     % 同义保留字段，便于脚本中显式说明

par.expClip = 700;
par.gaugeMode = 'max';                % 旧 split 路径：max-gauge；IF 路径可设 'auto'/'incremental'
par.postGaugeMode = 'max';
par.ifGaugeBalance = false;
par.ifGaugeStrategy = 'none';         % v13 推荐 'phase-peak'
par.ifRhsLogTarget = [];
par.ifRhsLogClip = 120;
par.ifSeedLogAbsMax = 600;
par.ifSeedLogRangeMax = 500;
par.rejectOnSeedLogRange = false;
par.ifRhsLogHardClip = 650;
par.clipNegativeW = true;
par.WFloor = 0;
par.monitorMatrixCondition = false;
par.enableSmallEpsMonitor = false;
par.monitorWarnThreshold = 80;

% H 的标量 gauge。旧版使用原始 H，因此默认 'none'。
% 小 epsilon IF 实验可设 'min'，去掉 H 的公共偏移。
par.hamiltonianGauge = 'none';

% theta 依赖 rephase 默认关闭。它保持 n 不变，但会改变 u_theta；只有在明确
% 排查 W 振幅指数跨度时才建议打开。
par.thetaRephase = false;
par.thetaRephaseStatistic = 'median';
par.thetaRephaseSmoothPasses = 0;
par.thetaRephaseMaxAbsShift = 0;
par.thetaRephaseUsePeakGauge = false;

par.ifMode = 'none';                  % v13 wkb-if-split 推荐 'h-only'
par.reactionScaleMaxAllowed = 1e8;
par.reactionScaleMinAllowed = 1e-8;
par.reactionScaleRetryLimit = 1e4;
par.adaptiveIfLogRange = false;
par.ifLogRangeMax = 80;
par.rhoBlowupThreshold = 1e100;
par.rhoLogRetryMax = Inf;             % 旧版默认不因 rho log 尺度拒步
par.enableStepRetry = false;
par.maxStepRetries = 12;
par.retryDtFactor = 0.5;
par.dtMinRetry = 1e-12;

%=========================================================================
% 6. 可选自适应时间步
%=========================================================================
par.adaptiveTimeStep = false;
par.adaptiveStrategy = 'auto';
par.adaptiveSafety = 0.5;
par.phaseCflSafety = 0.45;
par.dtMax = [];
par.dtMin = 0;

%=========================================================================
% 7. 停止、历史、快照、checkpoint 与进度输出
%=========================================================================
par.stopByResidual = true;
par.storeSnapshots = true;
par.snapshotTimes = [];

% 中间输出约定：
%   snapshot_<solver>_stepXXXX_tXX.mat 直接保存在 outdir 下；
%   checkpoint_<solver>_latest.mat 覆盖保存最近状态；
%   history_<solver>_latest.mat 保存轻量 history/diagnostics。
% 默认开启 checkpoint，但完整周期快照只在 snapshotTimes 或 snapshotEveryTime
% 有限时写出；这兼容旧脚本，也方便中途停止后恢复/后处理。
par.saveIntermediate = [];
par.snapshotEveryTime = Inf;
par.snapshotEverySteps = [];
par.saveLatestCheckpoint = true;
par.checkpointEveryTime = 0.5;
par.checkpointEverySteps = [];
par.saveHistoryCheckpoint = true;
par.snapshotSubdir = '';  % 保留字段兼容旧脚本；新版不再创建 snapshots 子文件夹。
par.checkpointFile = 'checkpoint_wkb_latest.mat';
par.historyCheckpointFile = 'history_wkb_latest.mat';
par.printIntermediateSaveMessage = false;
par.saveMatV73 = false;
par.saveFinalResult = true;
par.saveInitialAcceptedSnapshot = true;

par.verbose = true;
par.progressEveryPercent = 5;
par.progressEverySeconds = 10;
par.progressEverySteps = [];
par.historyEveryTime = 0.5;
par.historyEverySteps = [];

%=========================================================================
% 8. 在线图像监测控制
%=========================================================================
par.livePlot = false;
par.livePlotMode = 'monitor';
par.livePlotEveryTime = 0.5;
par.livePlotEverySteps = [];
par.livePlotSave = false;
par.livePlotPause = 0.0;

%=========================================================================
% 9. 系数函数
%=========================================================================
switch lower(profile)
    case {'baseline','legacy','paper'}
        theta_m = 0.35;
        K0 = 1.0;
        K1 = 0.5;
        Dmin = 0.2;
        D1 = 0.4;
        par.profile = 'baseline';
        par.modelVariant = 'legacy-minus-cos';
        par.theta_m = theta_m;
        par.K0 = K0;
        par.K1 = K1;
        par.Dmin = Dmin;
        par.D1 = D1;
        par.K_fun = @(x) K0 - K1 * cos(2*pi*x);  % 旧版/论文基准设置
        par.D_fun = @(theta) Dmin + D1 * (1 - cos(2*pi*(theta - theta_m)));

    case {'baseline-v13','v13-baseline','baseline_plus','baseline-plus'}
        theta_m = 0.35;
        K0 = 1.0;
        K1 = 0.5;
        Dmin = 0.2;
        D1 = 0.4;
        par.profile = 'baseline-v13';
        par.modelVariant = 'v13-plus-cos';
        par.theta_m = theta_m;
        par.K0 = K0;
        par.K1 = K1;
        par.Dmin = Dmin;
        par.D1 = D1;
        par.K_fun = @(x) K0 + K1 * cos(2*pi*x);  % 仅用于复现 v13 曾用设置
        par.D_fun = @(theta) Dmin + D1 * (1 - cos(2*pi*(theta - theta_m)));

    case 'heterogeneous'
        theta_m = 0.25;
        par.theta_m = theta_m;
        par.modelVariant = 'heterogeneous';
        par.K_fun = @(x) 1 + 0.35*cos(2*pi*x) + 0.20*cos(4*pi*x);
        par.D_fun = @(theta) 0.15 + 0.35*(1 - cos(2*pi*(theta - theta_m))) + ...
                             0.08*(1 - cos(4*pi*(theta - theta_m)));
    otherwise
        error('Unknown profile "%s".', profile);
end
end
