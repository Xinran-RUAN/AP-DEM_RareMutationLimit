function original_bio

%% �������    
Nth = 7; % theta������������
Nx = 7; % x����  
theta = linspace(0, 1, Nth +1); % �����theta
theta(end) = [];  
x = linspace(0, 1, Nx + 1); % �����x  
x = x';     
dtheta = theta(2) - theta(1);
dx = x(2) - x(1);
dt = 1e-4;   
T = 50; 
Ts = 1:50;
ks = 1;

%% ��ֵ׼��
[xx, thth] = meshgrid(x, [theta, 1]);
thinte = 0:0.001:1; 
[xxinte, ththinte] = meshgrid(x, thinte);

%% �������   
epsilon = 1e-2;     
D = 0.5 * sin(pi * theta - pi) + 1;  
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)  
plot(theta, D);     
%% ��ֵw(x,theta,t),u(theta,t) 
w = ones(Nx + 1, Nth); % theta �����ڱ߽�������theta in [0, 1-dtheta]
u = sin(2*pi * theta) - 1;
W = exp(1 .* (sin(pi .* theta) + sin(pi .* x)).^2);
u = 0 .* theta;     
ne = w .* exp(u/epsilon); % theta �����ڱ߽�������theta in [0, 1-dtheta]
t = 0; 

mkdir('original_data');        
path = './original_data/'; 

%% ĳЩ׼�����֣�������ѭ����
beta = epsilon^2 * dt / dtheta^2;
alpha = dt .* D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(Nx-1, 1), 0) .* (1+ 2 * alpha + 2 * beta) + ...
    - diag(ones(Nx-2, 1), 1) .* alpha + ...
    - diag(ones(Nx-2, 1), -1) .* alpha;
A_pre(1, 1, :) = 1 + alpha + 2 * beta;
A_pre(end, end, :) = 1 + alpha + 2 * beta;
A_sub = -beta .* diag(ones(Nx-1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(Nth-1, 1), 1), A_sub) + kron(diag(ones(Nth-1, 1), -1), A_sub);
A_con(end-Nx+2:end, 1:Nx-1) = A_sub;
A_con(1:Nx-1, end-Nx+2:end) = A_sub;
path = './original_data/';
%% ʱ���ݻ�
while t <= T    
         
    if abs(t-Ts(ks)) <10^(-7)
        save(strcat(path, 'ne_', num2str(epsilon), '_', num2str(t), '_', num2str(Nx), '_', num2str(Nth), '.mat'), 'ne');
        ks = ks + 1;
    end

    %% rho,���֣���ֵ���֣�����epsilon�ļ�С��
    %% �����ֵ���ֲ����û᲻�������⣬����Ҳ��ﲻ����rho��x�йأ���theta�޹�      
   rho = dtheta * sum(ne(:, 1:Nth), 2);
    %% ���Բ�ֵ
    % ne_aux = [ne, ne(:, 1)];
    % neinte = interp2(xx, thth, ne_aux', xxinte, ththinte, 'spline');
    % neinte = neinte';
    % rho = 0.001 * sum(neinte(:, 1:end-1), 2);
    
%     %% ��ͼ  
%     figure(1);   
%     plot(x, rho);
    % axis([0 1 0 25]);   
 
    %% solve the equation to obtain n(x_j,\theta_i,t_m+1)
    %%%ϵ������   
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    A = A_con - dt * kron(diag(ones(Nth, 1), 0), A_diag2);

    % �Ҷ���  
    b = ne; 
    b = reshape(b(2: Nx, :), [], 1);
    %��ǰEuler       
    ne_new = A \ b;               
    ne = reshape(ne_new, size(ne(2:Nx, :)));% ����
    ne = [ne(1, :); ne; ne(end, :)];
    t = t+dt;    
 
%     figure(2)
%     ntheta = dx * (0.5*ne(1, :) + sum(ne(2:Nx, :), 1) + 0.5*ne(end, :));
%     plot(theta, ntheta);
end        

end

