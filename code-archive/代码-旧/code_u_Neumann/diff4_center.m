function d4u = diff4_center(u, dx)
% diff4_center - 周期边界四阶中心差分求四阶导数
% u : 输入向量
% dx : 空间步长
% d4u : 输出四阶导数近似

N = length(u);
d4u = zeros(size(u));

for j = 1:N
    % 周期索引
    jp1 = mod(j, N) + 1;
    jp2 = mod(j+1, N) + 1;
    jm1 = mod(j-2, N) + 1;
    jm2 = mod(j-3, N) + 1;

    % 四阶中心差分公式
    d4u(j) = (u(jm2) - 4*u(jm1) + 6*u(j) - 4*u(jp1) + u(jp2)) / dx^4;
end

end


