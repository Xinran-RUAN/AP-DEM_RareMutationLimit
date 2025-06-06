 %% ������ɢ��
    L = 8;           % �ռ䷶Χ [-L, L]
    Nx = 32;         % �ռ��������  
    Ny = 32;
    Nz = 32;
    dx = 2*L/Nx;     % �ռ䲽��
    dy = 2*L/Ny;
    dz = 2*L/Nz;
    x = linspace(-L, L, Nx);  % �ռ�����
    y = linspace(-L, L, Ny);
    z = linspace(-L, L, Nz);
    [X, Y, Z] = meshgrid(x, y, z);  % ��ά����
    tau = 0.1;      % ʱ�䲽��
    
   %% ��������
    beta = 0;        % �����Բ���   
    lambda3 = 0.1;   % LHY��ϵ��

    % ���Ƴ�
    V = 0.5 * (X.^2 + Y.^2 + Z.^2);
    v = reshape(V, Nx*Ny*Nz, 1);

   %%  ��ʼ������
    Phi = (1/pi)^(3/2) * exp(-(X.^2 + Y.^2 + Z.^2)/ 2);   % ��̬   
    Phi = Phi / sqrt(sum(sum(sum(abs(Phi).^2) * dx) * dy )* dz); % ��һ��
    phi = reshape(Phi, Nx*Ny*Nz, 1);
    
    %% ��������
    max_iter = 1000; % ����������
    tol = 1e-8;      % �����ݲ�
    
    %% ���׵�����
    % ��x���׵��Ĳ�־����ʾ
    d_xx = spdiags([ones(Nx,1), -2*ones(Nx,1), ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
    D_xx = kron(kron(speye(Nz), speye(Ny)), d_xx);

    % ��y���׵��Ĳ�־����ʾ
    d_yy = spdiags([ones(Ny,1), -2*ones(Ny,1), ones(Ny,1)], [-1 0 1], Ny, Ny) / dy^2;
    D_yy = kron(kron(speye(Nz), d_yy), speye(Nx));

    % ��z���׵��Ĳ�־����ʾ
    d_zz = spdiags([ones(Nz,1), -2*ones(Nz,1), ones(Nz,1)], [-1 0 1], Nz, Nz) / dz^2;
    D_zz = kron(kron(d_zz, speye(Ny)), speye(Nx));

    %% Crank-Nicolson��CNFD��
    for iter = 1:max_iter
        % ������������|��|? + ��?|��|?��
        nonlinear_term = beta * abs(phi).^2 + lambda3 * abs(phi).^3;
        
        % ������ܶ�������-0.5?? + V + ��������
        H_matrix = 0.5*(D_xx + D_yy + D_zz) - spdiags(v + nonlinear_term, 0, Nx*Ny*Nz, Nx*Ny*Nz);
        
        % ����CNFD����ϵͳ������Ҷ���
        A = speye(Nx*Ny*Nz) - 0.5*tau*H_matrix;
        rhs = (speye(Nx*Ny*Nz) + 0.5*tau*H_matrix) * phi;
        
        % �������ϵͳ
        phi_new = A \ rhs;
        
        % ��һ��������
        norm_factor = sqrt(sum(abs(phi_new).^2) * dx * dy * dz);
        phi_new = phi_new / norm_factor;
        
        % ���������
        if norm(phi_new - phi, 'inf') < tol
            fprintf('�����ڵ� %d �ε���\n', iter);
            break;
        end
        phi = phi_new;
    end
    
    % ������
    if iter == max_iter
        fprintf('�ﵽ������������δ����\n');
    end
    
    Phi = reshape(phi, Nx, Ny, Nz);
    V = reshape(v, Nx, Ny, Nz);
    
    %% ���ƻ�̬������
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
    title('��̬������');
    grid on;
    colorbar;
    
    %% �����ݶ�
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

    %% �����뻯ѧ��
    E_kinetic = 0.5 * (abs(dphi_dx).^2 + abs(dphi_dy).^2 + abs(dphi_dz).^2);
    E_potential = V .* abs(Phi).^2;
    E_nonlinear = 0.5 * beta * abs(Phi).^4 + (lambda3 / 4) * abs(Phi).^5;
    E = sum(E_kinetic(:) + E_potential(:) + E_nonlinear(:)) * dx * dy * dz;

    mu_integrand = E_kinetic + E_potential + beta * abs(Phi).^4 + lambda3 * abs(Phi).^5;
    mu = sum(mu_integrand(:)) * dx * dy * dz;

    fprintf('����: %.6f\n', E);
    fprintf('��ѧ��: %.6f\n', mu);
    fprintf('���: %.10f\n', norm(phi_next - phi_current, 'inf'));

