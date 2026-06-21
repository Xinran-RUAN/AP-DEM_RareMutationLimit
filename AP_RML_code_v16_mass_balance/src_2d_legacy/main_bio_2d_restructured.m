function result = main_bio_2d_restructured(par)
%MAIN_BIO_2D_RESTRUCTURED 二维 FEM legacy WKB 求解器封装。
%   有限元矩阵和旧 WENO/FEM 核心仍保存在 src_2d_legacy 中；本函数
%   增加统一参数、进度输出和最终诊断，便于与 1D run 脚本保持一致。
if nargin < 1
    par = model.default_params_2d();
end
here = fileparts(mfilename('fullpath'));
old = pwd;
cleanup = onCleanup(@() cd(old));
cd(here);
addpath(genpath(here));

N_theta = par.Ntheta;
d_theta = 1/N_theta;
theta = 0:d_theta:1-d_theta;
prog = utils.progress_init('WKB-2D-legacy', par);

% 生成/读取单位方形区域上的三角网格。
mesh_triangle(200, par.h, 0, 1, 0, 1);
load mesh_divide_elements.txt;
load mesh_divide_nodes.txt;
mesh_divide_elements = mesh_divide_elements; %#ok<NODEF,ASGSL>
mesh_divide_nodes = mesh_divide_nodes; %#ok<NODEF,ASGSL>
element_num = size(mesh_divide_elements,1);

x = mesh_divide_nodes(:,1);
y = mesh_divide_nodes(:,2);
D = par.D_fun(theta);
K = par.K_fun(x,y);

% 初始相位峰值刻意远离 theta_m，用来检验长期选择过程。
u = -0.05*(1 - cos(2*pi*(theta - par.theta0)));
W = repmat(1 + 0.1*cos(2*pi*x).*cos(2*pi*y), 1, N_theta);
[W, u] = src.normalize_gauge(W, u, par.eps, 'max');

% 预处理有限元积分点和刚度矩阵。
[xi, eta, weight] = numerical_intergal(22);
B = prepare_part(par.eps, par.dt, d_theta, N_theta);
[K1, K2_x, K2_y, Con] = element_stiffness_prepare_2D(xi, eta, mesh_divide_elements, mesh_divide_nodes, element_num);
rho = solve_rho(u, W, theta, par.eps, d_theta, N_theta, 1);

t = 0;
step = 0;
residual = Inf;
history.t = [];
history.theta_wkb = [];
history.residual = [];
historyEverySteps = local_history_every_steps(par);
lastInfo = struct();

while t < par.T && step < par.maxSteps
    H = solve_H_fem(N_theta, mesh_divide_elements, mesh_divide_nodes, D, K, rho, xi, eta, weight, K1, K2_x, K2_y);
    uOld = u;
    WOld = W;
    u = solve_u(u, B, H, d_theta, par.dt);
    W = solve_w_fem(u, W, K, rho, D, H, par.eps, d_theta, par.dt, N_theta, ...
                    mesh_divide_elements, mesh_divide_nodes, xi, eta, weight, K1, K2_x, K2_y, Con);
    [W, u] = src.normalize_gauge(W, u, par.eps, 'max');
    rho = solve_rho(u, W, theta, par.eps, d_theta, N_theta, 1);
    t = t + par.dt;
    step = step + 1;
    residual = max(max(abs(u(:)-uOld(:)))/par.dt, max(abs(W(:)-WOld(:)))/(par.dt*(1+max(abs(W(:))))));

    if step == 1 || mod(step, historyEverySteps) == 0 || residual <= par.tol
        [~, iu] = max(u);
        history.t(end+1,1) = t;
        history.theta_wkb(end+1,1) = theta(iu);
        history.residual(end+1,1) = residual;
        lastInfo.theta_wkb = theta(iu);
        lastInfo.err_wkb_to_m = utils.periodic_distance(theta(iu), par.theta_m);
        lastInfo.rho_max = max(rho);
    end
    prog = utils.progress_update(prog, step, t, residual, lastInfo);

    if residual < par.tol
        break;
    end
end
[~, iu] = max(u);
diagnostics.theta_wkb = theta(iu);
diagnostics.err_wkb_to_m = utils.periodic_distance(theta(iu), par.theta_m);
diagnostics.rho_max = max(rho);
diagnostics.W_min = min(W(:));
result = struct('W',W,'u',u,'rho',rho,'D',D,'K',K,'theta',theta,'nodes',mesh_divide_nodes, ...
                'elements',mesh_divide_elements,'par',par,'history',history, ...
                'diagnostics',diagnostics,'t',t,'step',step,'residual',residual);
utils.progress_finish(prog, step, t, residual, diagnostics);
end

function historyEverySteps = local_history_every_steps(par)
if isfield(par, 'historyEverySteps') && ~isempty(par.historyEverySteps) && ...
        isnumeric(par.historyEverySteps) && isscalar(par.historyEverySteps) && par.historyEverySteps > 0
    historyEverySteps = max(1, round(par.historyEverySteps));
elseif isfield(par, 'historyEveryTime') && ~isempty(par.historyEveryTime) && ...
        isnumeric(par.historyEveryTime) && isscalar(par.historyEveryTime) && par.historyEveryTime > 0
    historyEverySteps = max(1, round(par.historyEveryTime / par.dt));
else
    historyEverySteps = max(1, round(0.1 / par.dt));
end
end
