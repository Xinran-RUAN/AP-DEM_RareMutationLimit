function original_bio_notime

%% �������  
Nth = 50; % theta������������
Nx = 20; % x���� 
theta = linspace(0, 1, Nth +1); % �����theta
theta(end) = [];  
x = linspace(0, 1, Nx + 1); % �����x  
x = x';     
dtheta = theta(2) - theta(1);
dx = x(2) - x(1);
 
%% ������� 
epsilon = 1e-1;     
D = 0.5 * sin(pi * theta - pi) + 1;  
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)  
plot(theta, D);     
%% ��ֵw(x,theta,t),u(theta,t) �Լ���ʼ�� H
ne = ones(Nx + 1, Nth); % theta �����ڱ߽�������theta in [0, 1-dtheta]
t = 0; 

%% ĳЩ׼�����֣�������ѭ����
beta = epsilon^2 / dtheta^2;
alpha = D ./ dx^2;
alpha = reshape(alpha, 1, 1, []);
A_pre = diag(ones(Nx-1, 1), 0) .* (2 * alpha + 2 * beta) + ...
    - diag(ones(Nx-2, 1), 1) .* alpha + ...
    - diag(ones(Nx-2, 1), -1) .* alpha;
A_pre(1, 1, :) = alpha + 2 * beta;
A_pre(end, end, :) = alpha + 2 * beta;
A_sub = -beta .* diag(ones(Nx-1, 1), 0);
A_cell = squeeze(num2cell(A_pre, [1, 2]));
A_con = blkdiag(A_cell{:});
A_con = A_con + kron(diag(ones(Nth-1, 1), 1), A_sub) + kron(diag(ones(Nth-1, 1), -1), A_sub);
A_con(end-Nx+2:end, 1:Nx-1) = A_sub;
A_con(1:Nx-1, end-Nx+2:end) = A_sub;
%% ʱ���ݻ�
tol = 1;
while tol > 1e-9  
    %% rho,���֣���ֵ���֣�����epsilon�ļ�С��
    %% �����ֵ���ֲ����û᲻�������⣬����Ҳ��ﲻ����rho��x�йأ���theta�޹�      
    rho = dtheta * sum(ne(:, 1:Nth), 2);
      
    %% ��ͼ  
    figure(1);   
    plot(x, rho);
    axis([0 1 0 25]);   

    %% solve the equation to obtain n(x_j,\theta_i,t_m+1)
    %%%ϵ������
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    A = A_con - kron(diag(ones(Nth, 1), 0), A_diag2);

    % �Ҷ���  
    b = reshape(zeros(Nx-1, Nth), [], 1);
    %��ǰEuler       
    ne_new = A \ b;  
    ne_new = reshape(ne_new, size(ne(2:Nx, :)));% ����
    ne_new = [ne_new(1, :); ne_new; ne_new(end, :)];
    
    tol = max(max(abs(ne-ne_new)))
    ne = ne_new;
end

end

