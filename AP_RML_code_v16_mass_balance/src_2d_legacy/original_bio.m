function original_bio
clear  
close all 
%% 网格参数       
Nth = 500; % theta方向的网格点数
Nx = 5; % x方向   
theta = linspace(0, 1, Nth +1); % 网格点theta
theta(end) = [];  
x = linspace(0, 1, Nx + 1); % 网格点x  
x = x';       
dtheta = theta(2) - theta(1);
dx = x(2) - x(1);
dt = 1e-3;   
T = 2; 
Ts = 0.1:0.1:2;
ks = 1;

marker = 1;

%% 问题参数   
eps = 1e-0;     
D = initial_D(0.5, theta);
K = 1 + 20 * (1 - 4 * (x - 0.5).^2).^8;  
figure(10)    
plot(theta, D);     
%% 初值w(x,theta,t),u(theta,t) 
[W, u] = initial_Wu(theta, x);
[W, u] = normalize_u(W, u, eps, theta);  
ne = W .* exp(u/eps); % theta 是周期边界条件，theta in [0, 1-dtheta]
t = 0; 

baseOutdir = utils.output_dir('legacy_original_bio');
outdir = utils.ensure_dir(fullfile(baseOutdir, sprintf('eps_%s', utils.time_token(eps))));

%% 某些准备部分，不放在循环里
A_con = prepare_part_ori(eps, dx, dt, dtheta, D, Nx, Nth);

%% 时间演化
while t <= T    
         
    if abs(t-Ts(ks)) <10^(-7)
        snapshotFile = utils.legacy_snapshot_file(outdir, 'legacy_direct', Nx, Nth, t);
        save(snapshotFile, 'ne', 'eps', 't', 'Nx', 'Nth', 'theta', 'x', 'marker');
        ks = ks + 1;
        
        %% 画图
        ntheta = dx * (0.5*ne(1, :) + sum(ne(2:Nx, :), 1) + 0.5*ne(end, :));
        plot(theta, ntheta);
    end  

    %% rho,积分，数值积分，随着epsilon的减小，
    %% 这个数值积分不晓得会不会有问题，精度也许达不到。rho与x有关，与theta无关      
    rho = solve_original_rho(ne, x, theta);
 
    %% solve the equation to obtain n(x_j,\theta_i,t_m+1)
    %%系数矩阵   
    A_diag2 = diag(K - rho, 0); 
    A = A_con - dt * kron(diag(ones(Nth, 1), 0), A_diag2);
 
    % 右端项   
    b = eps * ne; 
    b = reshape(b, [], 1);
    %向前Euler       
    ne_new = A \ b;               
    ne = reshape(ne_new, size(ne));% 更新
   
    t = t+dt;     
end             

summaryFile = utils.summary_file(outdir);
save(summaryFile, 'eps', 'Nx', 'Nth', 'T', 'dt', 'Ts', 'outdir');

end
