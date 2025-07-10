function phix_minus = weno5_left(phi, dx)
% WENO5 ����͸������ع� (���� Shu-Osher)
% --------------------------------------------
% phi        : 1D ���飬������������ȡֵ
% dx         : ������ (����)
% phix_minus : WENO5 ����͸���������ֵ (1D ����)

N = length(phi);              % �������
phix_minus = zeros(size(phi)); % ��ʼ����������
epsilon = 1e-0;                % ��ֹ����Ĳ���

for i = 1:N
    % ������������ (MATLAB 1-based)
    i_m3 = mod(i - 4, N) + 1;  % i-3
    i_m2 = mod(i - 3, N) + 1;  % i-2
    i_m1 = mod(i - 2, N) + 1;  % i-1
    i_p1 = mod(i, N) + 1;      % i+1
    i_p2 = mod(i + 1, N) + 1;  % i+2

    % ���
    a = (phi(i_m2) - phi(i_m3)) / dx;
    b = (phi(i_m1) - phi(i_m2)) / dx;
    c = (phi(i)    - phi(i_m1)) / dx;
    d = (phi(i_p1) - phi(i)) / dx;
    e = (phi(i_p2) - phi(i_p1)) / dx;

    % ������ģ�嵼������
    phi0 = (1/3)*a - (7/6)*b + (11/6)*c;
    phi1 = -(1/6)*b + (5/6)*c + (1/3)*d;
    phi2 = (1/3)*c + (5/6)*d - (1/6)*e;

    % �⻬��ָ�� ������󣡣�������������������������������������
    % ����������������������������������������������������a,b,c,d���������е�abcd,
    % �����������µĹ�ʽ
    IS0 = 13*(a - b)^2 + 3*(a - 3*b)^2;
    IS1 = 13*(b - c)^2 + 3*(b - c)^2;
    IS2 = 13*(c - d)^2 + 3*(3*c - d)^2;

    % ������Ȩ��
    alpha0 = 1 / (epsilon + IS0)^2;
    alpha1 = 6 / (epsilon + IS1)^2;
    alpha2 = 3 / (epsilon + IS2)^2;
    alpha_sum = alpha0 + alpha1 + alpha2;
    w0 = alpha0 / alpha_sum;
    w1 = alpha1 / alpha_sum;
    w2 = alpha2 / alpha_sum;

    % WENO ����͵���
    phix_minus(i) = w0 * phi0 + w1 * phi1 + w2 * phi2;
end
end
