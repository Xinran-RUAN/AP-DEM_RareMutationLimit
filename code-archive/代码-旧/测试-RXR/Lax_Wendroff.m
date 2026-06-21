% Lax-Wendroff
% 参数设置
a = 1;          % 波速
L = 10;         % 空间域长度
T = 6;          % 时间域长度
Nx = 100;       % 空间网格点数
Nt = 300;       % 时间网格点数
dx = L / (Nx - 1);  % 空间步长
dt = T / (Nt - 1);  % 时间步长

% 初始化网格
x = linspace(0, L, Nx)';
t = linspace(0, T, Nt);
u = zeros(Nx, Nt);

% 初始条件：高斯脉冲
u(:, 1) = exp(-10 * (x - L / 2).^2);

% Lax-Wendroff格式时间演化
for n = 1:Nt-1
    for i = 2:Nx-1
        u(i, n+1) = u(i, n) - (a * dt / (2 * dx)) * (u(i+1, n) - u(i-1, n)) ...
                    + (a^2 * dt^2 / (2 * dx^2)) * (u(i+1, n) - 2 * u(i, n) + u(i-1, n));
    end

    % 边界条件：周期性边界条件或零边界条件，这里使用周期性边界条件
    u(1, n+1) = u(Nx, n+1);
    u(Nx, n+1) = u(1, n+1);
end

% 可视化结果
for n = 1:Nt
    plot(x, u(:, n));
    xlim([0 L]);
    ylim([-1 1]);
    title(['Lax-Wendroff格式求解PDE： Time = ', num2str(t(n))]);
    pause(0.01);
end 

% 作者：北太天元卢朓 https://www.bilibili.com/read/cv33293897/?jump_opus=1 出处：bilibili