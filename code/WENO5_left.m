function pL = WENO5_left(u, dx)
    eps0 = 1e-6;
    N = length(u);
    u = [u(end-2:end); u; u(1:3)];
    pL = zeros(N, 1);
    for j = 1:N
        jj = j + 3;
        vmm = u(jj-2);
        vm  = u(jj-1);
        v   = u(jj);
        vp  = u(jj+1);
        vpp = u(jj+2);

        % 三个候选
        q0 = (1/3) * vmm - (7/6) * vm + (11/6) * v;
        q1 = -(1/6) * vm + (5/6) * v + (1/3) * vp;
        q2 = (1/3) * v + (5/6) * vp - (1/6) * vpp;

        % 平滑指标
        beta0 = (13/12)*(vmm - 2*vm + v)^2 + 0.25*(vmm - 4*vm + 3*v)^2;
        beta1 = (13/12)*(vm - 2*v + vp)^2 + 0.25*(vm - vp)^2;
        beta2 = (13/12)*(v - 2*vp + vpp)^2 + 0.25*(3*v - 4*vp + vpp)^2;

        a0 = 0.1 / (eps0 + beta0)^2;
        a1 = 0.6 / (eps0 + beta1)^2;
        a2 = 0.3 / (eps0 + beta2)^2;
        wsum = a0 + a1 + a2;
        w0 = a0 / wsum;
        w1 = a1 / wsum;
        w2 = a2 / wsum;

        pL(j) = (w0 * q0 + w1 * q1 + w2 * q2) / dx;
    end
end