function du = central_diff_symm(u, dx, order)
% 高阶中心差分，一阶导数（对称 stencil，周期边界）
% u     : 输入向量
% dx    : 空间步长
% order : 阶数（4 或 6）

N = length(u);
du = zeros(size(u));

if order == 4
    for j = 1:N
        jp1 = mod(j, N) + 1;
        jp2 = mod(j+1, N) + 1;
        jm1 = mod(j-2, N) + 1;
        jm2 = mod(j-3, N) + 1;
        du(j) = (-u(jp2) + 8*u(jp1) - 8*u(jm1) + u(jm2)) / (12 * dx);
    end

elseif order == 6
    for j = 1:N
        jp1 = mod(j, N) + 1;
        jp2 = mod(j+1, N) + 1;
        jp3 = mod(j+2, N) + 1;
        jm1 = mod(j-2, N) + 1;
        jm2 = mod(j-3, N) + 1;
        jm3 = mod(j-4, N) + 1;
        du(j) = (u(jm3) - 9*u(jm2) + 45*u(jm1) ...
                - 45*u(jp1) + 9*u(jp2) - u(jp3)) / (60 * dx);
    end

else
    error('Order must be 4 or 6');
end

end
