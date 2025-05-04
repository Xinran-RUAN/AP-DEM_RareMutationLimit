function main_bio
%% 网格参数   
Nth = 20; % theta方向网格步长1/Nth
Nx = 20; % x方向网格步长1/Nx
theta = linspace(0, 1, Nth +1); % 网格点theta
theta(end) = []; %周期边界，去掉theta(end)=1 
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x'; % theta为行向量，x为列向量
dtheta = theta(2) - theta(1); %theta步长
dx = x(2) - x(1); %空间步长
dt = 1e-2; %时间步长  
T = 50; % 终止时间
Ts = 1:50; % 存储Ts时刻的数据
ks = 1; % 与Ts有关
  
%% 插值准备
[xx, thth] = meshgrid(x, [theta, 1]);
thinte = 0:0.001:1; 
[xx2, ththinte] = meshgrid(x, thinte);
 
%% 问题参数  
eps = 1e-2;       
D = 0.5 * sin(pi * theta - pi) + 1;     
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)  
plot(theta, D);    

%% 初值w(x,theta,t),u(theta,t) 以及初始化 H(\theta)
w = ones(Nx+1, Nth); % theta \in [0, 1-dtheta]
u = sin(pi * theta) - 1;
H = zeros(1, Nth);
t = 0; 

%% 某些准备部分，不放在循环里
%% 一个theta对应一个D，对应1个特征值，so, 要么循环，要么三阶张量形式
% 待定吧
beta = eps^2 * dt / dtheta^2;
B = (1 + 2*beta) * diag(ones(Nth, 1), 0) +...
    - beta * diag(ones(Nth-1, 1), 1) +...
    - beta *  diag(ones(Nth-1, 1), -1);
B(1, Nth) = -beta;
B(Nth, 1) = -beta;
% 二阶导隐式，系数矩阵  
alpha = dt .* D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(Nx-1, 1), 0) .* (1+2*alpha+2*beta) + ...
    - diag(ones(Nx-2, 1), 1) .* alpha + ...
    - diag(ones(Nx-2, 1), -1) .* alpha; 
A_pre(1, 1, :) = 1+alpha+2 * beta;
A_pre(end, end, :) = 1+alpha+2 * beta;
A_sub = -beta .* diag(ones(Nx-1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(Nth-1, 1), 1), A_sub) + kron(diag(ones(Nth-1, 1), -1), A_sub);
A_con(end-Nx+2:end, 1:Nx-1) = A_sub;
A_con(1:Nx-1, end-Nx+2:end) = A_sub;
path = './data/';
%% 时间演化
while t <= T  
    
    if abs(t-Ts(ks)) <10^(-7)
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'), 'u');
        save(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'), 'w');
        ks = ks + 1;
    end
    %% rho,积分，数值积分，随着epsilon的减小，
    %% 这个数值积分不晓得会不会有问题，精度也许达不到。rho与x有关，与theta无关
    %%直接算       
 %   rho1 = dtheta * sum(w(:, 1:Nth) .* exp(u(1:Nth)/eps), 2);
    %% 插值
    u_aux = [u, u(:, 1)];
    w_aux = [w, w(:, 1)];
    uinte = interp1([theta, 1], u_aux, thinte, 'spline');
    winte = interp2(xx, thth, w_aux', xx2, ththinte, 'spline');
    winte = winte';   
    rho = (thinte(2)-thinte(1)) * sum(winte(:, 1:end-1) .* exp(uinte(1:end-1)/eps), 2);
      
    %% 画图   
    figure(1);    
    plot(x, rho);
    %% 画图  
    figure(2);  
    plot(theta, u);  
 %   axis([0 1 0 2]);   
 %   pause(0.001); 
    
    %% 特征值问题。系数矩阵C表示， D(theta),
    %% 一个theta对应一个D，对应1个特征值，so, 循环
    Tri_C = - 2 * diag(ones(Nx-1, 1), 0) +...
            + diag(ones(Nx-2, 1), 1) +...
            + diag(ones(Nx-2, 1), -1);  
    Tri_C(1, 1) = -1;  
    Tri_C(end, end) = -1;
    for kk = 1: Nth % 计算特征值
        C = - D(kk) ./ dx^2 .* Tri_C - diag(K(2:Nx) - rho(2: Nx));
        H(kk) = eigs(C, 1, 'smallestreal');
    end    
    
    figure(3)
    plot(theta, H);
    
    %% solve the second equation to obtain u(theta)，迎风格式？
    u_plus = [u(2:end), u(1)]; % 关于theta是周期边界
    u_minus = [u(end), u(1:end-1)];
    % 非线性项的正负处理：信息p^2从左来，即一阶导p>0,则左差分/后向差分；否则右差分
    pL = (u - u_minus) / dtheta; % 左差分 u_j - u_{j-1}
    pR = (u_plus - u) / dtheta; % 右差分 u_{j+1} - u_j
    pC = (u_plus - u_minus) / (2*dtheta);
    grad_sq = pL.^2 .* (pC >= 0) + pR.^2 .* (pC < 0);
    % 二阶导用隐式，一阶导用迎风,隐式Euler
    u_new = B\(u + dt * eps * (grad_sq - H))';  
    
    %% solve the first equation to obtain W_varepsilon (x, theta)
    % 计算 du/dtheta, 正负用来判断风向
    u = u_new'; % 更新  
    u_plus = [u(2:end), u(1)];
    u_minus = [u(end), u(1:end-1)];
    dudtheta = (u_plus - u_minus) / (2*dtheta);
  %
    w_theta_plus = w(:, [2:end, 1]); % 关于theta是周期边界条件
    w_theta_minus = w(:, [end, 1:end-1]);
    dw_dtheta_forward = (w_theta_plus - w) / dtheta; % 右差分，dudtheta<0，
    dw_dtheta_backward = (w - w_theta_minus) / dtheta; % 左差分, dudtheta>0
    dw_dtheta_upwind = dw_dtheta_backward .* (dudtheta >= 0) +...
                     + dw_dtheta_forward .* (dudtheta <0);
    
    coupling = 2 * eps * (dw_dtheta_upwind .* dudtheta);
    
    %%%系数矩阵
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    
    H_mat = reshape(H, [1, 1, Nth]);
    A_H_cell = squeeze(num2cell(H_mat.*eye(Nx-1), [1, 2]));
    A_H = blkdiag(A_H_cell{:});
    A = A_con - dt * kron(diag(ones(Nth, 1), 0), A_diag2) - dt * A_H;

    % 右端项  
    b = w + dt * (coupling); 
    b = reshape(b(2: Nx, :), [], 1);
    %向前Euler         
    w_new = A \ b;      
    w = reshape(w_new, size(w(2:Nx, :)));% 更新
    w = [w(1, :); w; w(end, :)];
    t = t+dt;  
    
end

end

