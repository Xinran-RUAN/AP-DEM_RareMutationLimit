% 测试如下方程的求解方法
%   u_t - |u_s|^2 - eps * u_ss = -H
% s \in [0,1], 周期边界
% H = (x - 0.4)^2
%% 参数设置
ds = 1e-2;
s = 0:ds:1-ds;
N = length(s);
H_func = @(s) (s- 0.4).^2;
H = H_func(s);

%% 矩阵构造
MAT_D = Matrix_Laplacian1D_Periodic(N, ds);
MAT_U = speye(N) + dt * eps * MAT_D;

%% 时间推进
t = 0;
while t < T
     if t + dt > T
        dt = T - t;
    end
    u = RK3_step(u, dx, dt, eps, H, MAT_U);
    t = t + dt;
end
