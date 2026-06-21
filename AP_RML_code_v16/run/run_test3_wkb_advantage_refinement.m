% RUN_TEST3_WKB_ADVANTAGE_REFINEMENT
% -------------------------------------------------------------------------
% 算例 3：WKB 粗 trait 网格 vs direct trait 加密。
% 唯一批量变量是 Ntheta；WKB 和 direct 可各自给列表。
% -------------------------------------------------------------------------

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
startup_AP_RML(false);

%% ========================= 用户接口：算例名称和求解器 =========================

cfg.caseName = 'test3_wkb_advantage_refinement';
cfg.solvers = {'wkb','direct'};

%% ========================= 用户接口：基本数值参数 =========================

cfg.profile = 'baseline';
cfg.eps = 1e-4;
cfg.Nx = 40;
cfg.Ntheta_wkb = 15;
cfg.Ntheta_direct = [16];   % 唯一批量变量：Ntheta
cfg.dt = 1e-4;
cfg.T = 2;
cfg.tol = 1e-10;

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

cfg.u0_fun = @(theta) - cfg.eps * (1 - cos(2*pi*(theta - 0.70))); % 初值 O(epsilon)
cfg.W0_fun = @(x, theta) ones(numel(x), numel(theta));

%% ========================= 用户接口：算法选择 =========================

cfg.timeIntegrator = 'frozen';              % {'frozen'}
cfg.phaseHamiltonian = 'weno5';             % {'weno5','godunov','lf'}
cfg.amplitudeVariant = 'split';             % {'split','density-compatible-semidiscrete','wkb-if-split','split-if'}
cfg.rhoReconstruction = 'direct-log';       % {'direct-log','laplace-hybrid','laplace'} 测试
cfg.reactionDiscretization = 'implicit';    % {'implicit','patankar','logistic-exact'}
cfg.residualMode = 'legacy';                % {'legacy','wkb-steady','scaled-gauge-invariant'}

cfg.transportImplicit = true;
cfg.hamiltonianGauge = 'none';
cfg.thetaRephase = false;
cfg.smallEpsAutoPatches = false;
cfg.enableV13AutoPatches = false;

%% ========================= 用户接口：保存和输出 =========================

cfg.snapshotTimes = unique([0, 0.5, 1, cfg.T]);
cfg.snapshotEveryTime = Inf;
cfg.saveFinalResult = true;
cfg.saveLatestCheckpoint = false;
cfg.saveHistoryCheckpoint = true;
cfg.historyEveryTime = 0.5;

cfg.clearOutput = true;
cfg.verbose = true;
cfg.progressEveryPercent = 5;
cfg.progressEverySeconds = 20;
cfg.livePlot = false;

%% ========================= 运行 =========================

runInfo = run_case_1d(cfg); %#ok<NASGU>


%% ========================= 算法选择 =========================
% 说明：
%   这里的选项只控制 WKB/direct 求解中的数值处理方式。
%   若当前核心程序暂时没有实现某个选项，运行时会在参数检查处报错。
%   一般复现实验或与旧版本对比时，建议先使用下面这一组默认组合。

% -------------------------------------------------------------------------
% 时间推进方式
% -------------------------------------------------------------------------
% 可选：
%   'frozen'
%       冻结部分系数的时间推进方式。
%       这是目前最稳妥、最接近旧版程序的选择，适合先复现实验。
%
%   'standard'
%       标准时间推进方式。
%       若核心程序中保留该选项，可用于和 frozen 策略作对比。
%
%   'projective'
%       投影型或加速型时间推进方式。
%       适合尝试长时间演化，但更依赖稳定性设置。
%


% -------------------------------------------------------------------------
% phase 方程中 Hamiltonian 或 theta 方向导数的离散方式
% -------------------------------------------------------------------------
% 可选：
%   'weno5'
%       五阶 WENO 重构。
%       适合处理 u(theta) 中可能出现的不光滑点或尖峰附近的非光滑结构。
%       用于后处理或高分辨率比较时，通常比 spline 更可靠。
%
%   'godunov'
%       Godunov 型数值 Hamiltonian。
%       单调性好，理论分析中常用，但精度相对低一些。
%
%   'lf'
%       Lax-Friedrichs 型数值 Hamiltonian。
%       更带数值黏性，通常较稳，但可能更耗散。
%
%   'central'
%       中心差分。
%       形式简单，但在 Hamilton-Jacobi 型方程中一般不作为主要推荐选项。
%


% -------------------------------------------------------------------------
% amplitude 变量 W 的离散格式
% -------------------------------------------------------------------------
% 可选：
%   'split'
%       分裂型 WKB amplitude 格式。
%       这是目前主要使用的版本，便于分别处理 phase 和 amplitude。
%
%   'density'
%       密度兼容型变体。
%       更强调重构密度 n = W exp(u/eps) 后的守恒或平衡结构。
%
%   'legacy'
%       旧版 amplitude 处理方式。
%       主要用于复现旧代码结果和排查新旧版本差异。
%


% -------------------------------------------------------------------------
% rho 或 trait marginal 的重构方式
% -------------------------------------------------------------------------
% 可选：
%   'direct'
%       直接用 theta 网格求和或梯形积分。
%       简单直接，但当 eps 很小时，trait 分布高度集中，粗网格可能不够准。
%
%   'quadrature'
%       基于离散网格的数值积分。
%       适合一般精度要求，也便于和理论中的离散求和对应。
%
%   'laplace'
%       Laplace 或 steepest-descent 型重构。
%       适合 rare mutation limit 下 eps 很小、分布强集中时的 rho 重构。
%       例 3 中若要体现 WKB 粗 theta 网格的优势，建议优先使用该选项。
%


% -------------------------------------------------------------------------
% reaction 项的离散方式
% -------------------------------------------------------------------------
% 可选：
%   'implicit'
%       隐式处理 reaction 项。
%       稳定性更好，尤其适合反应项较强或时间步不太小时。
%
%   'explicit'
%       显式处理 reaction 项。
%       实现简单，但稳定性对 dt 更敏感。
%
%   'semi-implicit'
%       半隐式处理。
%       若核心程序支持，可用于在稳定性和计算成本之间折中。
%


% -------------------------------------------------------------------------
% residual 或 remainder 的处理模式
% -------------------------------------------------------------------------
% 可选：
%   'legacy'
%       旧版残差处理方式。
%       用于复现旧程序结果，建议作为默认对照。
%
%   'ap'
%       面向 AP 结构的残差处理方式。
%       若核心程序已实现，可用于测试固定网格下 eps -> 0 的结构保持性质。
%
%   'none'
%       不额外处理 residual。
%       主要用于调试，不建议作为正式数值实验默认选项。
%

