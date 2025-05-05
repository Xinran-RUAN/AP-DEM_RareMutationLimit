%% �������   
N_theta = 40; % theta�������񲽳�1/Nth
N_x = 20; % x�������񲽳�1/Nx
d_theta = 1 / N_theta;
theta = 0:d_theta:1-d_theta;
dx = 1 / N_x;
x = (0:dx:1)'; % �����x  
dt = 1e-1; %ʱ�䲽��  
T = 10000; % ��ֹʱ��
Ts = 1:50; % �洢Tsʱ�̵�����
ks = 1; % ��Ts�й�
  
%% ��ֵ׼��
[X, Theta] = meshgrid(x, [theta, 1]);
theta_f = 0:0.001:1; 
[X_f, Theta_f] = meshgrid(x, theta_f);
 
%% �������  
eps = 1e-5;       
D = 0.5 * sin(pi * theta - pi) + 1;     
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  

%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H(\theta)
W = ones(N_x+1, N_theta); % theta \in [0, 1-dtheta]
u = 100*(sin(2*pi * theta) - 1);
H = zeros(1, N_theta);
t = 0;    

%% ׼�����֣�������ѭ����
%% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, Ҫôѭ����Ҫô����������ʽ
% ������
beta = eps^2 * dt / d_theta^2;  
beta2 = eps^2 * dt / d_theta^2;
B = (1 + 2*beta2) * diag(ones(N_theta, 1), 0) +...
    - beta2 * diag(ones(N_theta-1, 1), 1) +...
    - beta2 *  diag(ones(N_theta-1, 1), -1);
B(1, N_theta) = -beta2;
B(N_theta, 1) = -beta2;
% ���׵���ʽ��ϵ������  
alpha = dt .* D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(N_x-1, 1), 0) .* (1+2*alpha+2*beta) + ...
    - diag(ones(N_x-2, 1), 1) .* alpha + ...
    - diag(ones(N_x-2, 1), -1) .* alpha; 
A_pre(1, 1, :) = 1+alpha+2 * beta;
A_pre(end, end, :) = 1+alpha+2 * beta;
A_sub = -beta .* diag(ones(N_x-1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(N_theta-1, 1), 1), A_sub) + kron(diag(ones(N_theta-1, 1), -1), A_sub);
A_con(end-N_x+2:end, 1:N_x-1) = A_sub;
A_con(1:N_x-1, end-N_x+2:end) = A_sub;
path = './data/';
%% ʱ���ݻ�
while t <= T  
    
    if min(abs(t-Ts)) < dt/10
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '.mat'), 'u');
        save(strcat(path, 'W_', num2str(eps), '_', num2str(t), '_', num2str(N_x), '_', num2str(N_theta), '.mat'), 'W');
    end
    %% rho,���֣���ֵ���֣�����epsilon�ļ�С��
    %% �����ֵ���ֲ����û᲻�������⣬����Ҳ��ﲻ����rho��x�йأ���theta�޹�
    %%ֱ����       
    rho = d_theta * sum(W(:, 1:N_theta) .* exp(u(1:N_theta)/eps), 2);
    %% ��ֵ
%     u_aux = [u, u(:, 1)];
%     w_aux = [W, W(:, 1)];
%     uinte = interp1([theta, 1], u_aux, theta_f, 'spline');
%     winte = interp2(X, Theta, w_aux', X_f, Theta_f, 'spline');
%     winte = winte';   
%     rho = (theta_f(2)-theta_f(1)) * sum(winte(:, 1:end-1) .* exp(uinte(1:end-1)/eps), 2);
%       
    %% ��ͼ   
    figure(1);    
    plot(x, rho);
    title(['Plot of $\rho(x)$, time = ', num2str(t)], 'Interpreter', 'latex');
    %% ��ͼ  
    figure(2);    
    plot(theta, u); 
    title(['Plot of $u(\theta)$, time = ', num2str(t)], 'Interpreter', 'latex');
    
    %% ����ֵ���⡣ϵ������C��ʾ�� D(theta),
    %% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, ѭ��
    Tri_C = - 2 * diag(ones(N_x-1, 1), 0) +...
            + diag(ones(N_x-2, 1), 1) +...
            + diag(ones(N_x-2, 1), -1);  
    Tri_C(1, 1) = -1;  
    Tri_C(end, end) = -1;
    for kk = 1: N_theta % ��������ֵ
        C = - D(kk) ./ dx^2 .* Tri_C - diag(K(2:N_x) - rho(2: N_x));
        H(kk) = eigs(C, 1, 'smallestreal');
    end    
    
    figure(3)
    plot(theta, H);
    title(['Plot of $H(\theta)$, time = ', num2str(t)], 'Interpreter', 'latex');

    %% solve the second equation to obtain u(theta)��ӭ���ʽ��
    u_plus = [u(2:end), u(1)]; % ����theta�����ڱ߽�
    u_minus = [u(end), u(1:end-1)];
    % �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
    pL = (u - u_minus) / d_theta; % ���� u_j - u_{j-1}
    pR = (u_plus - u) / d_theta; % �Ҳ�� u_{j+1} - u_j
    pC = (u_plus - u_minus) / (2*d_theta);
    grad_sq = pL.^2 .* (pC >= 0) + pR.^2 .* (pC < 0);
    grad_sq = (pL .* (pC>=0) + pR.*(pC<0)) .* pC;
    % ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
    u_new = B\(u + dt * eps * (grad_sq - H))';  
    
    %% solve the first equation to obtain W_varepsilon (x, theta)
    % ���� du/dtheta, ���������жϷ���
    u = u_new'; % ����    
    u_plus = [u(2:end), u(1)];
    u_minus = [u(end), u(1:end-1)];
    dudtheta = (u_plus - u_minus) / (2*d_theta);
  %
    w_theta_plus = W(:, [2:end, 1]); % ����theta�����ڱ߽�����
    w_theta_minus = W(:, [end, 1:end-1]);
    dw_dtheta_forward = (w_theta_plus - W) / d_theta; % �Ҳ�֣�dudtheta<0��
    dw_dtheta_backward = (W - w_theta_minus) / d_theta; % ����, dudtheta>0
    dw_dtheta_upwind = dw_dtheta_backward .* (dudtheta >= 0) +...
                     + dw_dtheta_forward .* (dudtheta <0);
    
    coupling = 2 * eps * (dw_dtheta_upwind .* dudtheta);
    
    %%%ϵ������
    A_diag2 = diag(K(2:N_x) - rho(2:N_x), 0); 
    
    H_mat = reshape(H, [1, 1, N_theta]);
    A_H_cell = squeeze(num2cell(H_mat.*eye(N_x-1), [1, 2]));
    A_H = blkdiag(A_H_cell{:});
    A = A_con - dt * kron(diag(ones(N_theta, 1), 0), A_diag2) - dt * A_H;

    % �Ҷ���  
    b = W + dt * (coupling); 
    b = reshape(b(2: N_x, :), [], 1);
    %��ǰEuler         
    w_new = A \ b;      
    W = reshape(w_new, size(W(2:N_x, :)));% ����
    W = [W(1, :); W; W(end, :)];
    t = t+dt;  
    
end  
