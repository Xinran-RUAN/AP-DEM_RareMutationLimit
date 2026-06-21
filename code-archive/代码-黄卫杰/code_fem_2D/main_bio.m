clear  
close all                  
       
%% ol = 1代表spline插值；否则代表未spline插值
% spline插值的地方：u的方程中u对theta的一阶导；
% W方程中，u对theta的一阶导，W对theta的一阶导
ol = 0;            
        
%% 网格参数                     
N_theta = 15; % theta方向网格步长1/Nth
d_theta = 1 / N_theta;      
theta = 0:d_theta:1-d_theta;% 周期边界条件

%% 时间参数
tau = 1e-3; %时间步长 
T = 2; % 终止时间         
Ts = 0.1:0.1:T; % 存储数据的时刻  
  
%% 空间三角形网格剖分
%% 网格剖分，均匀剖分
h = 0.1; 
xa = 0; 
xb = 1;  
yc = 0; 
yd = 1;
mesh_triangle(200, h, xa, xb, yc, yd);
%% 下载数据  
load mesh_divide_elements.txt;
load mesh_divide_boundary_nodes.txt;
load mesh_divide_nodes.txt;
%% 网格剖分信息
element_num = length(mesh_divide_elements);
node_num = length(mesh_divide_nodes);
%% 判断哪些节点位于边界
nodes_boundary = find(mesh_divide_boundary_nodes == 1);

%% 问题参数             
eps = 1e-3;        
D = initial_D(0.5, theta); 
plot(theta, D);  
x = mesh_divide_nodes(:, 1);
y = mesh_divide_nodes(:, 2);
K = 1 + 20 * (1 - 2 * (x - 0.5).^2 - 2 * (y - 0.5).^2).^8;  
mkdir('data');            
path = './data/';   
plot_K(x, y, K); 
   
%% 初值w(x,theta,t),u(theta,t) 以及初始化 H(\theta)
[W, u] = initial_Wu(theta, x, y);
plot_K(x, y, W(:, 5));
[W, u] = normalize_u(W, u, eps, theta);
t = 0;        
tol = 10^(-12);
[xi, eta, weight] = numerical_intergal(22);  

%% 准备部分，为了节省时间，不放在循环里
B = prepare_part(eps, tau, d_theta, N_theta);
[K1, K2_x, K2_y, Con] = element_stiffness_prepare_2D(xi, eta, mesh_divide_elements, mesh_divide_nodes, element_num);
rho = solve_rho(u, W, theta, eps, d_theta, N_theta, ol);
%% 时间演化    
while t <= T       
    %% 存储数据 
    dt = tau;           
    if min(abs(t-Ts)) < dt/100
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '.mat'), 'u');
        save(strcat(path, 'W_', num2str(eps), '_', num2str(t), '.mat'), 'W');
          plot_all(mesh_divide_elements, mesh_divide_nodes, theta, rho, t, u, W, H, eps, 0);
    end      

    
    %% 特征值问题。  
    H = solve_H_fem(N_theta, mesh_divide_elements,mesh_divide_nodes,...
                    D, K, rho, xi, eta, weight, K1, K2_x, K2_y);
              
    %% solve the second equation to obtain u(theta)
    u_new = solve_u(u, B, H, d_theta, dt);  
    u = u_new;    

    W_new = solve_w_fem(u, W, K, rho, D, H, eps, d_theta, dt, N_theta,...
                        mesh_divide_elements, mesh_divide_nodes,...
                        xi, eta, weight, K1, K2_x, K2_y, Con);   
    W = W_new;
    [W, u] = normalize_u(W, u, eps, theta);
  
    t = t+dt;    
    
   %% solve rho   
    rho = solve_rho_laplace(eps, theta, u, W, d_theta);  
   %% 画图          
 
 %  plot_all(mesh_divide_elements, mesh_divide_nodes, theta, rho, t, u, W, H, eps, 0)
end    