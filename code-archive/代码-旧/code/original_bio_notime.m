function original_bio_notime

%% 网格参数  
Nth = 50; % theta方向的网格点数
Nx = 20; % x方向 
theta = linspace(0, 1, Nth +1); % 网格点theta
theta(end) = [];  
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';     
dtheta = theta(2) - theta(1);
dx = x(2) - x(1);
 
%% 问题参数 
epsilon = 1e-1;     
D = 0.5 * sin(pi * theta - pi) + 1;  
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)  
plot(theta, D);     
%% 初值w(x,theta,t),u(theta,t) 以及初始化 H
ne = ones(Nx + 1, Nth); % theta 是周期边界条件，theta in [0, 1-dtheta]
t = 0; 

%% 某些准备部分，不放在循环里
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
%% 时间演化
tol = 1;
while tol > 1e-9  
    %% rho,积分，数值积分，随着epsilon的减小，
    %% 这个数值积分不晓得会不会有问题，精度也许达不到。rho与x有关，与theta无关      
    rho = dtheta * sum(ne(:, 1:Nth), 2);
      
    %% 画图  
    figure(1);   
    plot(x, rho);
    axis([0 1 0 25]);   

    %% solve the equation to obtain n(x_j,\theta_i,t_m+1)
    %%%系数矩阵
    A_diag2 = diag(K(2:Nx) - rho(2:Nx), 0); 
    A = A_con - kron(diag(ones(Nth, 1), 0), A_diag2);

    % 右端项  
    b = reshape(zeros(Nx-1, Nth), [], 1);
    %向前Euler       
    ne_new = A \ b;  
    ne_new = reshape(ne_new, size(ne(2:Nx, :)));% 更新
    ne_new = [ne_new(1, :); ne_new; ne_new(end, :)];
    
    tol = max(max(abs(ne-ne_new)))
    ne = ne_new;
end

end

