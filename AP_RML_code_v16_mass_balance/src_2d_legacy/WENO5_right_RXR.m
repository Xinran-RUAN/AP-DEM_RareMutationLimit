function pR = WENO5_right_RXR(u, dx)
u = u';
    eps0 = 1e-6;
    N = length(u);
    u = [u(end-2:end); u; u(1:3)];

    pR = zeros(N,1);
    for j = 1:N
        jj = j + 3;
        vppp = u(jj+3);
        vpp = u(jj+2);
        vp  = u(jj+1);
        v   = u(jj);
        vm  = u(jj-1);
        vmm = u(jj-2);

        %% 以下为点值的重构，�?�非单侧导数的重�?
        % 三个候�?? 
        % q0(j) = -(2*vmm - 7*vm + 11*v)/6;
        % q1(j) = -(-vm + 5*v + 2*vp)/6;
        % q2(j) = -(2*v + 5*vp - vpp)/6;


        % 平滑指标
        % beta0 = (13/12)*(vmm - 2*vm + v)^2 + 0.25*(vmm - 4*vm + 3*v)^2;
        % beta1 = (13/12)*(vm - 2*v + vp)^2 + 0.25*(vm - vp)^2;
        % beta2 = (13/12)*(v - 2*vp + vpp)^2 + 0.25*(3*v - 4*vp + vpp)^2;
        % 
        % a0 = 0.1 / (eps0 + beta0)^2;
        % a1 = 0.6 / (eps0 + beta1)^2;
        % a2 = 0.3 / (eps0 + beta2)^2;
        % wsum = a0 + a1 + a2;
        % w0 = a0 / wsum;
        % w1 = a1 / wsum;
        % w2 = a2 / wsum;

        % pR(j) = (w0 * q0(j) + w1 * q1(j) + w2 * q2(j)) / dx;
        %% 单侧导数的重�?
        a = vppp - 2 * vpp + vp;
        b = vpp - 2 * vp + v;
        c = vp - 2 * v + vm;
        d = v - 2 * vm + vmm;

        IS0 = 13 * (a - b)^2 + 3 * (a - 3 * b)^2;
        IS1 = 13 * (b - c)^2 + 3 * (b + c)^2;
        IS2 = 13 * (c - d)^2 + 3 * (3 * c - d)^2;

        a0 = 1 / (eps0 + IS0)^2;
        a1 = 6 / (eps0 + IS1)^2;
        a2 = 3 / (eps0 + IS2)^2;
        wsum = a0 + a1 + a2;
        w0 = a0 / wsum;
        % w1 = a1 / wsum;
        w2 = a2 / wsum;

        pR(j) = (...
            - (vm - vmm) ...
            + 7 * (v - vm) ...
            + 7 * (vp - v) ...
            - (vpp - vp)...
            ) / 12 ...
            + w0 * (a - 2 * b + c) / 3 + (w2 - 1/2) * (b - 2 * c + d) / 6;
        pR(j) = pR(j) / dx;
    end
    pR = pR';
end
