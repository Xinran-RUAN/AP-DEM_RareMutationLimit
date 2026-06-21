function pR = WENO5_right(u, dx)
    u = reshape(u, [length(u), 1]);
    % 周期边界
    u = [u(end-2:end); u; u(1:3)];
    
    % 取 stencil
    N = length(u) - 6;
    pR = zeros(N, 1);
    
    for j = 1:N
        vmm = u(j+4);
        vm  = u(j+3);
        v   = u(j+2);
        vp  = u(j+1);
        vpp = u(j);
        
        % 三个候选
        q1 = -(2*vmm - 7*vm + 11*v)/6;
        q2 = -(-vm + 5*v + 2*vp)/6;
        q3 = -(2*v + 5*vp - vpp)/6;
        
        % 平滑指标
        beta1 = 13/12*(vmm - 2*vm + v)^2 + 1/4*(vmm - 4*vm + 3*v)^2;
        beta2 = 13/12*(vm - 2*v + vp)^2 + 1/4*(vm - vp)^2;
        beta3 = 13/12*(v - 2*vp + vpp)^2 + 1/4*(3*v - 4*vp + vpp)^2;
        
        % 权重
        eps_weno = 1e-6;
        alpha1 = 0.1 ./ (eps_weno + beta1).^2;
        alpha2 = 0.6 ./ (eps_weno + beta2).^2;
        alpha3 = 0.3 ./ (eps_weno + beta3).^2;
        alpha_sum = alpha1 + alpha2 + alpha3;
        
        w1 = alpha1 / alpha_sum;
        w2 = alpha2 / alpha_sum;
        w3 = alpha3 / alpha_sum;
        
        % 重构
        pR(j) = (w1 * q1 + w2 * q2 + w3 * q3) / dx;
    end
end
