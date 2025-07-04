clear  
close all        
  
%% ol = 1����spline��ֵ���������δspline��ֵ
% spline��ֵ�ĵط���u�ķ�����u��theta��һ�׵���
% W�����У�u��theta��һ�׵���W��theta��һ�׵�
ol = 0;       
marker = 0;      
%% �������                 
N_theta = 20; % theta�������񲽳�1/Nth
N_x = 10; % x�������񲽳�1/Nx  
 
d_theta = 1 / N_theta;      
theta = 0:d_theta:1-d_theta;

d_x = 1 / N_x;  
x = (0:d_x:1)'; % ����� x       
 
%% ʱ�����
tau = 1e-3; %ʱ�䲽�� 
T = 50; % ��ֹʱ��        
Ts = 0.1:0.1:50; % �洢���ݵ�ʱ��  

%% �������             
eps = 1e-2;      
D = initial_D(0.6, theta);
plot(theta, D);   
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
mkdir('data');           
path = './data/';       
    
%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H(\theta)
[W, u] = initial_Wu(theta, x);
[W, u] = normalize_u(W, u, eps, theta);
t = 0;    
tol = 10^(-12);

%% ׼�����֣�Ϊ�˽�ʡʱ�䣬������ѭ����
[B, Ap, Tri_C] = prepare_part(eps, tau, d_x, d_theta, N_x, N_theta, D);
% B Ϊ���u��ϵ�������׼����AΪ���W��ϵ�������׼��
MAT_D = Matrix_Laplacian1D_Periodic(N_theta, d_theta);
MAT_U = speye(N_theta) + tau * eps * MAT_D;

%% ʱ���ݻ�    
while t <= T     
    %% �洢���� 
    dt = tau;           
    if min(abs(t-Ts)) < dt/100
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '_', num2str(ol), '_', num2str(marker), '.mat'), 'u');
        save(strcat(path, 'W_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '_', num2str(ol), '_', num2str(marker), '.mat'), 'W');
        plot_all(x, theta, rho, t, d_x, u, W, H, eps, N_x, 0)
    end      
    %% solve rho  
    rho = solve_rho(u, W, x, theta, eps, d_theta, N_theta, marker);
      
    %% ����ֵ���⡣
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
   %% ��ͼ     

     % plot_all(x, theta, rho, t, d_x, u, W, H, eps, N_x, 0)
end   