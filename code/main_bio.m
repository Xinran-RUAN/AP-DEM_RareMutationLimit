clear  
close all        
  
%% ol = 1代表spline插值；否则代表未spline插值
% spline插值的地方：u的方程中u对theta的一阶导；
% W方程中，u对theta的一阶导，W对theta的一阶导
ol = 0;       
marker = 0;      
%% 网格参数                 
N_theta = 20; % theta方向网格步长1/Nth
N_x = 10; % x方向网格步长1/Nx  
 
d_theta = 1 / N_theta;      
theta = 0:d_theta:1-d_theta;

d_x = 1 / N_x;  
x = (0:d_x:1)'; % 网格点 x       
 
%% 时间参数
tau = 1e-3; %时间步长 
T = 50; % 终止时间        
Ts = 0.1:0.1:50; % 存储数据的时刻  

%% 问题参数             
eps = 1e-2;      
D = initial_D(0.6, theta);
plot(theta, D);   
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
mkdir('data');           
path = './data/';       
    
%% 初值w(x,theta,t),u(theta,t) 以及初始化 H(\theta)
[W, u] = initial_Wu(theta, x);
[W, u] = normalize_u(W, u, eps, theta);
t = 0;    
tol = 10^(-12);

%% 准备部分，为了节省时间，不放在循环里
[B, Ap, Tri_C] = prepare_part(eps, tau, d_x, d_theta, N_x, N_theta, D);
% B 为求解u的系数矩阵的准备，A为求解W的系数矩阵的准备
MAT_D = Matrix_Laplacian1D_Periodic(N_theta, d_theta);
MAT_U = speye(N_theta) + tau * eps * MAT_D;

%% 时间演化    
while t <= T     
    %% 存储数据 
    dt = tau;           
    if min(abs(t-Ts)) < dt/100
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '_', num2str(ol), '_', num2str(marker), '.mat'), 'u');
        save(strcat(path, 'W_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '_', num2str(ol), '_', num2str(marker), '.mat'), 'W');
        plot_all(x, theta, rho, t, d_x, u, W, H, eps, N_x, 0)
    end      
    %% solve rho  
    rho = solve_rho(u, W, x, theta, eps, d_theta, N_theta, marker);
      
    %% 特征值问题。
    H = solve_H(N_theta, N_x, D, Tri_C, K, rho);
       
    %% solve the second equation to obtain u(theta)
     % u_new = solve_u(u, u, B, H, d_theta, dt, theta, ol); 
   u_new = solve_u_RXR(u, MAT_U, H, d_theta, dt); 
    
    if norm(u-u_new, 2) < tol
        break;
    end      
    u = u_new;    
    W = solve_w(u, W, W, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta, theta, ol);  
         
    [W, u] = normalize_u(W, u, eps, theta);


    t = t+dt;   
    if abs(t-0.4) < 1e-7
        u;
    end  
   %% 画图     

     % plot_all(x, theta, rho, t, d_x, u, W, H, eps, N_x, 0)
end   