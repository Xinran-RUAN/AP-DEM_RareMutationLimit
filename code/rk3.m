function phi = rk3(phi, dx, dt)
    k1 = weno5_hj_rhs(phi, dx);
    k2 = weno5_hj_rhs(phi + dt * k1 / 2, dx);
    k3 = weno5_hj_rhs(phi - dt * k1 + 2 * dt * k2, dx);
    
    phi = phi + dt/6 * (k1 + 4*k2 + k3);
end
