function du = central_diff_symm(u, dx, order)
% �߽����Ĳ�֣�һ�׵������Գ� stencil�����ڱ߽磩
% u     : ��������
% dx    : �ռ䲽��
% order : ������4 �� 6��

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
