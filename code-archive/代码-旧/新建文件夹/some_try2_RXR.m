%% 网格参数   
N_theta = 20; % theta方向网格步长1/Nth
N_x = 10; % x方向网格步长1/Nx
d_theta = 1 / N_theta;  
theta = 0:d_theta:1-d_theta;      
dx = 1 / N_x;        
x = (0:dx:1)'; % 网格点x  
dt = 1e-3; %时间步长       
T = 1000; % 终止时间  
Ts = 1:30; % 存储Ts时刻的数据
ks = 1; % 与Ts有关 
    
%% 插值准备
[X, Theta] = meshgrid(x, [theta, 1]);
theta_f = 0:0.001:1; 
[X_f, Theta_f] = meshgrid(x, theta_f);
   
%% 问题参数      
eps = 1e-8;           
D = 0.5 * sin(pi * theta - pi) + 1;    
D = exp(-1 .* sin(pi .* theta).^2); 
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  

%% 初值w(x,theta,t),u(theta,t) 以及初始化 H(\theta)
u = exp(1 .* (sin(pi * theta) + 0. * sin(2*pi*theta)).^2);
u = u- max(u);  
t = 0; 

%% 一个theta对应一个D，对应1个特征值，so, 要么循环，要么三阶张量形式
% 待定吧      
beta = eps * dt / d_theta^2;        
B = (1 + 2 * beta) * diag(ones(N_theta, 1), 0) +...
            - beta * diag(ones(N_theta-1, 1), 1) +...
            - beta * diag(ones(N_theta-1, 1), -1);
B(1, N_theta) = - beta;
B(N_theta, 1) = - beta;  
path = './data/';   
%% 时间演化  
while t <= T  
    
    %% rho,积分，数值积分，随着epsilon的减小，
    %% 这个数值积分不晓得会不会有问题，精度也许达不到。rho与x有关，与theta无关
  %% 画图    
    figure(2);     
    plot([theta, 1], [u, u(1)]);   
    title(['Plot of $u(\theta)$, time = ', num2str(t)], 'Interpreter', 'latex');

%     figure(3);  
%     plot(theta(1:end-1), diff(u));   
%     title(['Plot of $diff(u)$, time = ', num2str(t)], 'Interpreter', 'latex');
%       
    %% solve the second equation to obtain u(theta)，迎风格式？
    u_plus = [u(2:end), u(1)]; % 关于theta是周期边界
    u_minus = [u(end), u(1:end-1)];
    % 非线性项的正负处理：信息p^2从左来，即一阶导p>0,则左差分/后向差分；否则右差分
    pL = (u - u_minus) / d_theta; % 左差分 u_j - u_{j-1}
    pR = (u_plus - u) / d_theta; % 右差分 u_{j+1} - u_j
    grad_sq = min(pL.^2, pR.^2) .* (pL .* pR >= 0) + ...
              pL .* 0 .* (pL .* pR < 0); 
    grad_sq(1) = pR(1).^2;
    
%% 
    % 二阶导用隐式，一阶导用迎风,隐式Euler
    u_new = B\(u + dt * grad_sq)';  
    u = u_new';   
    %% 显示Euler呢？  
%     laplace_u = (u_plus - 2 * u + u_minus) / d_theta^2;
%     u_new = u + dt * (grad_sq + eps * laplace_u);
%     u = u_new;
    t=t+dt;     
    
end   
