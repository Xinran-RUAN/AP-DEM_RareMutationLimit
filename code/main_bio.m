function main_bio
%% �������   
Nth = 20; % theta�������񲽳�1/Nth
Nx = 20; % x�������񲽳�1/Nx
theta = linspace(0, 1, Nth +1); % �����theta
theta(end) = []; %���ڱ߽磬ȥ��theta(end)=1 
x = linspace(0, 1, Nx + 1); % �����x  
x = x'; % thetaΪ��������xΪ������
dtheta = theta(2) - theta(1); %theta����
dx = x(2) - x(1); %�ռ䲽��
dt = 1e-2; %ʱ�䲽��  
T = 50; % ��ֹʱ��
Ts = 1:50; % �洢Tsʱ�̵�����
ks = 1; % ��Ts�й�
  
%% ��ֵ׼��
[xx, thth] = meshgrid(x, [theta, 1]);
thinte = 0:0.001:1; 
[xx2, ththinte] = meshgrid(x, thinte);
 
%% �������  
eps = 1e-2;       
D = 0.5 * sin(pi * theta - pi) + 1;     
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)  
plot(theta, D);    

%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H(\theta)
w = ones(Nx+1, Nth); % theta \in [0, 1-dtheta]
u = sin(pi * theta) - 1;
H = zeros(1, Nth);
t = 0; 

%% ĳЩ׼�����֣�������ѭ����
%% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, Ҫôѭ����Ҫô����������ʽ
% ������
beta = eps^2 * dt / dtheta^2;
B = (1 + 2*beta) * diag(ones(Nth, 1), 0) +...
    - beta * diag(ones(Nth-1, 1), 1) +...
    - beta *  diag(ones(Nth-1, 1), -1);
B(1, Nth) = -beta;
B(Nth, 1) = -beta;
% ���׵���ʽ��ϵ������  
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
%% ʱ���ݻ�
while t <= T  
    
    if abs(t-Ts(ks)) <10^(-7)
        save(strcat(path, 'u_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'), 'u');
        save(strcat(path, 'w_', num2str(eps), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'), 'w');
        ks = ks + 1;
    end
    %% rho,���֣���ֵ���֣�����epsilon�ļ�С��
    %% �����ֵ���ֲ����û᲻�������⣬����Ҳ��ﲻ����rho��x�йأ���theta�޹�
    %%ֱ����       
 %   rho1 = dtheta * sum(w(:, 1:Nth) .* exp(u(1:Nth)/eps), 2);
    %% ��ֵ
    u_aux = [u, u(:, 1)];
    w_aux = [w, w(:, 1)];
    uinte = interp1([theta, 1], u_aux, thinte, 'spline');
    winte = interp2(xx, thth, w_aux', xx2, ththinte, 'spline');
    winte = winte';   
    rho = (thinte(2)-thinte(1)) * sum(winte(:, 1:end-1) .* exp(uinte(1:end-1)/eps), 2);
      
    %% ��ͼ   
    figure(1);    
    plot(x, rho);
    %% ��ͼ  
    figure(2);  
    plot(theta, u);  
 %   axis([0 1 0 2]);   
 %   pause(0.001); 
    
    %% ����ֵ���⡣ϵ������C��ʾ�� D(theta),
    %% һ��theta��Ӧһ��D����Ӧ1������ֵ��so, ѭ��
    Tri_C = - 2 * diag(ones(Nx-1, 1), 0) +...
            + diag(ones(Nx-2, 1), 1) +...
            + diag(ones(Nx-2, 1), -1);  
    Tri_C(1, 1) = -1;  
    Tri_C(end, end) = -1;
    for kk = 1: Nth % ��������ֵ
        C = - D(kk) ./ dx^2 .* Tri_C - diag(K(2:Nx) - rho(2: Nx));
        H(kk) = eigs(C, 1, 'smallestreal');
    end    
    
    figure(3)
    plot(theta, H);
    
    %% solve the second equation to obtain u(theta)��ӭ���ʽ��
    u_plus = [u(2:end), u(1)]; % ����theta�����ڱ߽�
    u_minus = [u(end), u(1:end-1)];
    % �������������������Ϣp^2����������һ�׵�p>0,������/�����֣������Ҳ��
    pL = (u - u_minus) / dtheta; % ���� u_j - u_{j-1}
    pR = (u_plus - u) / dtheta; % �Ҳ�� u_{j+1} - u_j
    pC = (u_plus - u_minus) / (2*dtheta);
    grad_sq = pL.^2 .* (pC >= 0) + pR.^2 .* (pC < 0);
    % ���׵�����ʽ��һ�׵���ӭ��,��ʽEuler
    u_new = B\(u + dt * eps * (grad_sq - H))';  
    
    %% solve the first equation to obtain W_varepsilon (x, theta)
    % ���� du/dtheta, ���������жϷ���
    u = u_new'; % ����  
    u_plus = [u(2:end), u(1)];
    u_minus = [u(end), u(1:end-1)];
    dudtheta = (u_plus - u_minus) / (2*dtheta);
  %
    w_theta_plus = w(:, [2:end, 1]); % ����theta�����ڱ߽�����
    w_theta_minus = w(:, [end, 1:end-1]);
    dw_dtheta_forward = (w_theta_plus - w) / dtheta; % �Ҳ�֣�dudtheta<0��
    dw_dtheta_backward = (w - w_theta_minus) / dtheta; % ����, dudtheta>0
    dw_dtheta_upwind = dw_dtheta_backward .* (dudtheta >= 0) +...
                     + dw_dtheta_forward .* (dudtheta <0);
    
    coupling = 2 * eps * (dw_dtheta_upwind .* dudtheta);
    
    %%%ϵ������
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    
    H_mat = reshape(H, [1, 1, Nth]);
    A_H_cell = squeeze(num2cell(H_mat.*eye(Nx-1), [1, 2]));
    A_H = blkdiag(A_H_cell{:});
    A = A_con - dt * kron(diag(ones(Nth, 1), 0), A_diag2) - dt * A_H;

    % �Ҷ���  
    b = w + dt * (coupling); 
    b = reshape(b(2: Nx, :), [], 1);
    %��ǰEuler         
    w_new = A \ b;      
    w = reshape(w_new, size(w(2:Nx, :)));% ����
    w = [w(1, :); w; w(end, :)];
    t = t+dt;  
    
end

end

