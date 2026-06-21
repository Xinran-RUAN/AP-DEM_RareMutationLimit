% 这里每一步都解隐式的线性系统，同时显式处理一阶项。
function u_new = RK3_step(u, dx, dt, H, MAT_U)
    rhs1 = RHS(u, dx, H);
    u1 = MAT_U \ (u + dt * rhs1);

    rhs2 = RHS(u1, dx, H);
    u2 = MAT_U \ ((3*u + u1 + dt * rhs2) / 4);

    rhs3 = RHS(u2, dx, H);
    u_new = MAT_U \ ((u + 2*u2 + 2*dt * rhs3) / 3);
end
