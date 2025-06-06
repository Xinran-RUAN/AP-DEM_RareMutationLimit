%% �������   
N_theta = 20; % theta�������񲽳�1/Nth
N_x = 10; % x�������񲽳�1/Nx
d_theta = 1 / N_theta;  
theta = 0:d_theta:1-d_theta;      
dx = 1 / N_x;        
x = (0:dx:1)'; % �����x  
dt = 1e-3; %ʱ�䲽��       
T = 1000; % ��ֹʱ��  
Ts = 1:30; % �洢Tsʱ�̵�����
ks = 1; % ��Ts�й� 
    
%% ��ֵ׼��
[X, Theta] = meshgrid(x, [theta, 1]);
theta_f = 0:0.001:1; 
[X_f, Theta_f] = meshgrid(x, theta_f);
   
%% �������      
eps = 1e-8;           
D = 0.5 * sin(pi * theta - pi) + 1;    
D = exp(-1 .* sin(pi .* theta).^2); 
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  

%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H(\theta)
u = exp(1 .* (sin(pi * theta) + 0. * sin(2*pi*theta)).^2);
u = u- max(u);  
t = 0; 

%% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, Ҫôѭ����Ҫô����������ʽ
% ������      
beta = eps * dt / d_theta^2;        
B = (1 + 2 * beta) * diag(ones(N_theta, 1), 0) +...
            - beta * diag(ones(N_theta-1, 1), 1) +...
            - beta * diag(ones(N_theta-1, 1), -1);
B(1, N_theta) = - beta;
B(N_theta, 1) = - beta;  
path = './data/';   
%% ʱ���ݻ�  
while t <= T  
    
    %% rho,���֣���ֵ���֣�����epsilon�ļ�С��
    %% �����ֵ���ֲ����û᲻�������⣬����Ҳ��ﲻ����rho��x�йأ���theta�޹�
  %% ��ͼ    
    figure(2);     
    plot([theta, 1], [u, u(1)]);   
    title(['Plot of $u(\theta)$, time = ', num2str(t)], 'Interpreter', 'latex');

%     figure(3);  
%     plot(theta(1:end-1), diff(u));   
%     title(['Plot of $diff(u)$, time = ', num2str(t)], 'Interpreter', 'latex');
%       
    %% solve the second equation to obtain u(theta)��ӭ���ʽ��
    u_plus = [u(2:end), u(1)]; % ����theta�����ڱ߽�
    u_minus = [u(end), u(1:end-1)];
    % �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
    pL = (u - u_minus) / d_theta; % ���� u_j - u_{j-1}
    pR = (u_plus - u) / d_theta; % �Ҳ�� u_{j+1} - u_j
    grad_sq = min(pL.^2, pR.^2) .* (pL .* pR >= 0) + ...
              pL .* 0 .* (pL .* pR < 0); 
    grad_sq(1) = pR(1).^2;
    
%% 
    % ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
    u_new = B\(u + dt * grad_sq)';  
    u = u_new';   
    %% ��ʾEuler�أ�  
%     laplace_u = (u_plus - 2 * u + u_minus) / d_theta^2;
%     u_new = u + dt * (grad_sq + eps * laplace_u);
%     u = u_new;
    t=t+dt;     
    
end   
