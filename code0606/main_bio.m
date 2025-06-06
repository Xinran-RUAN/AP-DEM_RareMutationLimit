clear  
close all  
%% �������         
N_theta = 7; % theta�������񲽳�1/Nth
N_x = 7; % x�������񲽳�1/Nx
d_theta = 1 / N_theta;      
theta = 0:d_theta:1-d_theta;
dx = 1 / N_x;
x = (0:dx:1)'; % �����x     
         
dt0 = 1e-4; %ʱ�䲽��              
T = 10; % ��ֹʱ��      
Ts = 1:1:T; % �洢Tsʱ�̵����� 
Con = 2;

%% �������           
eps = 1e-2;                   
%D = 0.5 * sin(pi * theta - pi) + 0.6; 
D = exp(-1 .* sin(pi .* theta).^2); % theta_m = 0.5;
D(2:4) = D(end:-1:end-2);
% D = exp(-1 .* (sin(pi * theta) + 0.5 * sin(2 * pi * theta)).^2); % theta_m \neq 0.5
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
if ~exist('data', 'dir')
    mkdir('data');
end
path = './data/'; 
  
%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H(\theta)
W = ones(N_x+1, N_theta); % theta \in [0, 1-dtheta]
% u = (sin(pi * theta) - 1);
u = exp(1 .* (sin(pi * theta) + 0. * sin(2*pi*theta)).^2);
W = exp(1 .* (sin(pi .* theta) + sin(pi .* x)).^2);
u = 0 .* theta;     
u = u - max(u);             
t = 0;  
tol = 10^(-12);
  
%% ׼�����֣�Ϊ�˽�ʡʱ�䣬������ѭ����
%% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, Ҫôѭ����Ҫô����������ʽ
%% B Ϊ���u��ϵ�������׼����AΪ���W��ϵ�������׼��

%% ʱ���ݻ�    
% for kk = 1: 1000
   [B, Ap, Tri_C] = prepare_part(eps, dt0, dx, d_theta, N_x, N_theta, D);
%    t = 0;
while t <= T    
    %% �洢���� 
    dt = dt0;       
    if min(abs(t-Ts)) < dt/100
        % save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '.mat'), 'u');
        % save(strcat(path, 'W_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '.mat'), 'W');
        plot_all(x, theta, rho, t, dx, u, W, H, eps, N_x, 0)
    end 
    %% solve rho
    rho = solve_rho(u, W, x, theta, eps, d_theta, N_theta);
    
    %% ����ֵ���⡣ϵ������C��ʾ�� D(theta),
    %% һ��theta��Ӧһ��D(theta)����Ӧ1������ֵ��������Ҫѭ��
    H = solve_H(N_theta, N_x, D, Tri_C, K, rho);
      
    %% solve the second equation to obtain u(theta)
    u_new = solve_u(u, u, B, H, d_theta, dt); 
    if norm(u-u_new, 2) < tol
        break;
    end     
    u = u_new;
    W = solve_w(u, W, W, K, rho, H, Ap, eps, d_theta, dt, N_x, N_theta);
 
 %   err = V ./ W(2:end-1, :)
       
    t = t+dt;
    % if min(abs(t-1:T)) < dt/5        
    %     disp(t)
    % end
   %% ��ͼ     
   % plot_all(x, theta, rho, t, dx, u, W, H, eps, N_x, 0)
    
end   
% kk = kk + 1;
% eps = eps * 0.5;
% end
