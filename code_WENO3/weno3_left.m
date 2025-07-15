function phix_minus = weno3_left(phi, dx)
% WENO3 ����͸������ع�
N = length(phi);
phix_minus = zeros(size(phi));
epsilon = 1e-1;

for i = 1:N
    i_m2 = mod(i-3, N) + 1;
    i_m1 = mod(i-2, N) + 1;
    i_p0 = i;
    i_p1 = mod(i, N) + 1;

    % ���
    a = (phi(i_p0) - phi(i_m1)) / dx; % ��_i - ��_{i-1}
    b = (phi(i_m1) - phi(i_m2)) / dx; % ��_{i-1} - ��_{i-2}
    c = (phi(i_p1) - phi(i_p0)) / dx; % ��_{i+1} - ��_i

    % ��ģ�嵼��
    phi0 = (3/2)*a - (1/2)*b;  % stencil {i-2,i-1,i}
    phi1 = (1/2)*a + (1/2)*c;  % stencil {i-1,i,i+1}

    % �⻬��ָ��
    IS0 = (phi(i_p0) - 2*phi(i_m1) + phi(i_m2))^2;
    IS1 = (phi(i_p1) - 2*phi(i_p0) + phi(i_m1))^2;

    % Ȩ��
    alpha0 = (1/3) / (epsilon + IS0)^2;
    alpha1 = (2/3) / (epsilon + IS1)^2;
    alpha_sum = alpha0 + alpha1;
    w0 = alpha0 / alpha_sum;
    w1 = alpha1 / alpha_sum;

    % �ϳɵ���
    phix_minus(i) = w0 * phi0 + w1 * phi1;
end
end
