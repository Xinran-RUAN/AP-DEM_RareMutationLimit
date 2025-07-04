clear  
close all
%% �������         
N_theta = 20; % theta�������񲽳�1/Nth
N_x = 20; % x�������񲽳�1/Nx
d_theta = 1 / N_theta;   
theta = 0:d_theta:1-d_theta;
dx = 1 / N_x;
x = (0:dx:1)'; % �����x    
      
dt0 = 1e-9; %ʱ�䲽��            
T = 10000; % ��ֹʱ��        
Ts = 1:50; % �洢Tsʱ�̵�����
tol = 10^(-10);   

%% �������           
eps = 1e-10;   
Con = 2;

%D = 0.5 * sin(pi * theta - pi) + 0.6; 
D = exp(-1 .* sin(pi .* theta).^2); % theta_m = 0.5;
% D = exp(-1 .* (sin(pi * theta) + 0.5 * sin(2 * pi * theta)).^2); % theta_m \neq 0.5
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
mkdir('data');        
path = './data/'; 
  
%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H(\theta)
W= ones(N_x + 1, N_theta);
% W = 0.1*exp(1 .* (sin(pi .* theta) + sin(pi .* x)).^2);
% W = exp(-sin(pi .* theta).^2 - sin(pi.*x).^2); % ones(N_x+1, N_theta) ; % theta \in [0, 1-dtheta]
% W = W - max(W);
% u = (sin(pi * theta) - 1);   
u = exp(1 .* (sin(pi * theta) + 0. * sin(2*pi*theta)).^2);
u = u - max(u);     
W = exp(1 .* (sin(pi .* theta) + sin(pi .* x)).^2);
u = 0 .* theta;
t = 0;          
    
%% ׼�����֣�Ϊ�˽�ʡʱ�䣬������ѭ����
%% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, Ҫôѭ����Ҫô����������ʽ
%% B Ϊ���u��ϵ�������׼����AΪ���W��ϵ�������׼��

%% ʱ���ݻ�  
while t <= T   
    dt = dt0;
    [B, Ap, Tri_C] = prepare_part(eps, dt, dx, d_theta, N_x, N_theta, D);

    %% �洢����
    if min(abs(t-Ts)) < dt/100
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '.mat'), 'u');
        save(strcat(path, 'W_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '.mat'), 'W');
    end 
    err = 1;
%     ne = W .* exp(u/eps);
%     ne = solve_n(eps, dt, dx, d_theta, N_x, N_theta, ne, D, K);
%     rho = d_theta * sum(ne(:, 1:N_theta), 2);
    rho = solve_rho(u, W, x, theta, eps);
    u_temp = u;
    W_temp = W;    
    while err > tol  
        %% solve rho   
         [~, H] = solve_H(N_theta, N_x, D, Tri_C, K, rho);
         u_temp = solve_u(u, u_temp, B, H, d_theta, dt); 
           
%          max_u = max(u_temp2);
%          while max_u > Con * eps
%              dt = dt / 2; 
%              [B, Ap, Tri_C] = prepare_part(eps, dt, dx, d_theta, N_x, N_theta, D);
%              u_temp2 = solve_u(u, u_temp, B, H, d_theta, dt); 
%              max_u = max(u_temp2);
%          end 
%          u_temp = u_temp2;
         % u_temp(u_temp>0) = 0;
         % u_temp = u_temp - max(u_temp);  
         W_temp = solve_w(u_temp, W, W_temp, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta);
         rho_temp = solve_rho(u_temp, W_temp, x, theta, eps); 
           
         err = norm(rho - rho_temp)
         rho = rho_temp; 
             
          plot_all(x, theta, rho_temp, t, dx, u_temp, W_temp, H, eps, N_x, 0)
    end 
    u = u_temp; 
    W = W_temp;
    t = t+dt;    
    
    %% ��ͼ 
    plot_all(x, theta, rho, t, dx, u, W, H, eps, N_x)
end  
