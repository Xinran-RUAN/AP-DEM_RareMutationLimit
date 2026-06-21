 %% 网格离散化
    L = 16;           % 空间范围 [-L, L]
    Nx = 32;         % 空间网格点数
    Ny = 32;   
    Nz = 32;   
    dx = 2*L/Nx;     % 空间步长
    dy = 2*L/Ny;
    dz = 2*L/Nz;  
    x = linspace(-L, L, Nx);  % 空间网格
    y = linspace(-L, L, Ny);
    z = linspace(-L, L, Nz);
    [X, Y, Z] = meshgrid(x, y, z);  % 三维网格
    tau = 0.1;      % 时间步长
    
   %% 参数设置  
    beta = 0;        % 非线性参数   
    lambda3 = 0.1;   % LHY项系数

    % 外势场
    V = 0.5 * (X.^2 + Y.^2 + Z.^2);
    v = reshape(V, Nx*Ny*Nz, 1);  

   %%  初始波函数
    % Phi = (1/pi)^(3/2) * exp(-(X.^2 + Y.^2 + Z.^2)/ 2);   % 基态   
    Phi = max(1 - X.^2 - Y.^2 - Z.^2, 0);   
    Phi = Phi / sqrt(sum(sum(sum(abs(Phi).^2) * dx) * dy )* dz); % 归一化
    phi_prev = reshape(Phi, Nx*Ny*Nz, 1);  % 初始 n-1
    phi_current = phi_prev;                % 初始 n
      
    %% 迭代参数
    max_iter = 1000;  % 最大迭代次数
    tol = 1e-10;      % 收敛                                                                                                                                                                     容差
    
    %% 二阶导矩阵
    % 对x二阶导的差分矩阵表示
    d_xx = spdiags([ones(Nx,1), -2*ones(Nx,1), ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
    D_xx = kron(kron(speye(Nz), speye(Ny)), d_xx);

    % 对y二阶导的差分矩阵表示
    d_yy = spdiags([ones(Ny,1), -2*ones(Ny,1), ones(Ny,1)], [-1 0 1], Ny, Ny) / dy^2;
    D_yy = kron(kron(speye(Nz), d_yy), speye(Nx));

    % 对z二阶导的差分矩阵表示
    d_zz = spdiags([ones(Nz,1), -2*ones(Nz,1), ones(Nz,1)], [-1 0 1], Nz, Nz) / dz^2;
    D_zz = kron(kron(d_zz, speye(Ny)), speye(Nx));

    %% SIFD
    % 初始处理
    nonlinear_term = beta * abs(phi_current).^2 + lambda3 * abs(phi_current).^3;
    H_linear = 0.5*(D_xx + D_yy + D_zz) - spdiags(v, 0, Nx*Ny*Nz, Nx*Ny*Nz);
    % 中间
    phi_predict = phi_current - 0.5 * tau * (H_linear * phi_current + nonlinear_term .* phi_current);
    % 第一步
    nonlinear_term = beta * abs(phi_predict).^2 + lambda3 * abs(phi_predict).^3;
    phi_next = phi_current - tau * (H_linear* phi_predict + nonlinear_term.*phi_predict);
    phi_prev = phi_current;
    phi_current = phi_next / sqrt(sum(abs(phi_next).^2) * dx * dy * dz); % 归一化
    
    %% 矩阵求解循环
    for iter = 1:max_iter
        % 计算非线性项
        nonlinear_term = beta * abs(phi_current).^2 + lambda3 * abs(phi_current).^3;
        
        % 构造线性系统矩阵
        A = speye(Nx*Ny*Nz) - tau * H_linear;
        rhs = (speye(Nx*Ny*Nz) + tau * H_linear) * phi_prev - 2*tau * nonlinear_term .* phi_current;
        % 求解线性系统
        phi_next = A \ rhs;   
               
        % 归一化  
        norm_factor = sqrt(sum(abs(phi_next).^2) * dx * dy * dz);
        phi_next = phi_next / norm_factor;
        
        % 检查收敛性
        if norm(phi_next - phi_current, 'inf') < tol
            fprintf('收敛于第 %d 次迭代\n', iter);
            break;
        end
        
        % 更新
        phi_prev = phi_current;
        phi_current = phi_next;
        
    Phi = reshape(phi_current, Nx, Ny, Nz);
    E = compute_energy(Phi, V, Nx, Ny, Nz, dx, dy, dz, beta, lambda3);
    figure(1);
    plot(x, Phi(:,Ny/2+1,Nz/2+1), 'LineWidth', 1);
    title(strcat('iter=', num2str(iter), ',   energy = ', num2str(E)));
    xlabel('x');
    ylabel('phi');
    grid on;
    axis([-L L -0.1 0.5]);
    end
    
    % 输出结果
    if iter == max_iter
        fprintf('达到最大迭代次数，未收敛\n');
    end
    
    Phi = reshape(phi_current, Nx, Ny, Nz);
    V = reshape(v, Nx, Ny, Nz);
    
    %% 绘图与能量计算
    figure(1);
    plot(x, Phi(:,Ny/2+1,Nz/2+1), 'LineWidth', 2);
    xlabel('x');
    ylabel('phi');
    grid on;
    figure(2);
    surf(x, y, Phi(:,:,Nz/2+1), 'EdgeColor', 'none');
    xlabel('x');
    ylabel('y');
    zlabel('phi');
    title('基态波函数');
    grid on;
    colorbar;
    
     %% 计算梯度
    dphi_dx = zeros(Nx,Ny,Nz);
    dphi_dx(2:Nx-1,:,:)=Phi(3:Nx,:,:)-Phi(1:Nx-2,:,:)/(2*dx);
    dphi_dx(1,:,:)=Phi(2,:,:)-Phi(1,:,:)/(dx);
    dphi_dx(Nx,:,:)=Phi(Nx,:,:)-Phi(Nx-1,:,:)/(dx);
    
    dphi_dy = zeros(Nx,Ny,Nz);
    dphi_dy(2:Ny-1,:,:)=Phi(3:Ny,:,:)-Phi(1:Ny-2,:,:)/(2*dy);
    dphi_dy(1,:,:)=Phi(2,:,:)-Phi(1,:,:)/(dy);
    dphi_dy(Nx,:,:)=Phi(Ny,:,:)-Phi(Ny-1,:,:)/(dy);
    
    dphi_dz = zeros(Nx,Ny,Nz);
    dphi_dz(2:Nz-1,:,:)=Phi(3:Nz,:,:)-Phi(1:Nz-2,:,:)/(2*dz);
    dphi_dz(1,:,:)=Phi(2,:,:)-Phi(1,:,:)/(dz);
    dphi_dz(Nx,:,:)=Phi(Nz,:,:)-Phi(Nz-1,:,:)/(dz);

   %% 能量与化学势
    E_kinetic = 0.5 * (abs(dphi_dx).^2 + abs(dphi_dy).^2 + abs(dphi_dz).^2);
    E_potential = V .* abs(Phi).^2;
    E_nonlinear = 0.5 * beta * abs(Phi).^4 + (lambda3 / 4) * abs(Phi).^5;
    E = sum(E_kinetic(:) + E_potential(:) + E_nonlinear(:)) * dx * dy * dz;

    mu_integrand = E_kinetic + E_potential + beta * abs(Phi).^4 + lambda3 * abs(Phi).^5;
    mu = sum(mu_integrand(:)) * dx * dy * dz;

    fprintf('能量: %.10f\n', E);
    fprintf('化学势: %.6f\n', mu);
    fprintf('误差: %.15f\n', norm(phi_next - phi_current, 'inf'));
