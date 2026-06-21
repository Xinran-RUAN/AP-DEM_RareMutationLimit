function u_new = solve_u_RXR(u, B, H, d_theta, dt)
    u = reshape(u, [length(u), 1]);
    H = reshape(H, [length(H), 1]);
    rhs = RHS(u, d_theta, H);
    u_new = B \ (u + dt * rhs);
    u_new = u_new';
end

