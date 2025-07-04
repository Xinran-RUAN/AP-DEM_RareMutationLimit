function original_bio
clear 
close all
%% �������       
Nth = 1000; % theta������������
Nx = 10; % x����   
theta = linspace(0, 1, Nth +1); % �����theta
theta(end) = [];  
x = linspace(0, 1, Nx + 1); % �����x  
x = x';       
dtheta = theta(2) - theta(1);
dx = x(2) - x(1);
dt = 1e-3;   
T = 10; 
Ts = 0.1:0.1:10;
ks = 1;

marker = 1;

%% �������   
eps = 1e-2;     
%D = 0.5 * sin(pi * theta - pi) + 1; 
D = initial_D(0.6, theta);
% D = exp(-1 .* (sin(pi * theta) + 0.4 * sin(2 * pi * theta)).^2); % theta_m \neq 0.5
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)   
plot(theta, D);     
%% ��ֵw(x,theta,t),u(theta,t) 
[W, u] = initial_Wu(theta, x);
[W, u] = normalize_u(W, u, eps, theta);  
ne = W .* exp(u/eps); % theta �����ڱ߽�������theta in [0, 1-dtheta]
t = 0; 

mkdir('original_data_nomid');        
path = './original_data_nomid/'; 

%% ĳЩ׼�����֣�������ѭ����
beta = eps^2 * dt / dtheta^2;
alpha = dt .* D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(Nx-1, 1), 0) .* (eps+ 2 * alpha + 2 * beta) + ...
    - diag(ones(Nx-2, 1), 1) .* alpha + ...
    - diag(ones(Nx-2, 1), -1) .* alpha;
A_pre(1, 1, :) = eps + alpha + 2 * beta;
A_pre(end, end, :) = eps + alpha + 2 * beta;
A_sub = -beta .* diag(ones(Nx-1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(Nth-1, 1), 1), A_sub) + kron(diag(ones(Nth-1, 1), -1), A_sub);
A_con(end-Nx+2:end, 1:Nx-1) = A_sub;
A_con(1:Nx-1, end-Nx+2:end) = A_sub;

%% ʱ���ݻ�
while t <= T    
         
    if abs(t-Ts(ks)) <10^(-7)
        save(strcat(path, 'ne_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '_', num2str(marker), '.mat'), 'ne');
        ks = ks + 1;
        
        %% ��ͼ
        ntheta = dx * (0.5*ne(1, :) + sum(ne(2:Nx, :), 1) + 0.5*ne(end, :));
        plot(theta, ntheta);
    end

    %% rho,���֣���ֵ���֣�����epsilon�ļ�С��
    %% �����ֵ���ֲ����û᲻�������⣬����Ҳ��ﲻ����rho��x�йأ���theta�޹�      
    rho = solve_original_rho(ne, Nth, x, theta);
 
    %% solve the equation to obtain n(x_j,\theta_i,t_m+1)
    %%ϵ������   
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    A = A_con - dt * kron(diag(ones(Nth, 1), 0), A_diag2);

    % �Ҷ���   
    b = eps * ne; 
    b = reshape(b(2: Nx, :), [], 1);
    %��ǰEuler       
    ne_new = A \ b;               
    ne = reshape(ne_new, size(ne(2:Nx, :)));% ����
    ne = [ne(1, :); ne; ne(end, :)];
    t = t+dt;    
end           

end
