function d4u = diff4_center(u, dx)
% diff4_center - ���ڱ߽��Ľ����Ĳ�����Ľ׵���
% u : ��������
% dx : �ռ䲽��
% d4u : ����Ľ׵�������

N = length(u);
d4u = zeros(size(u));

for j = 1:N
    % ��������
    jp1 = mod(j, N) + 1;
    jp2 = mod(j+1, N) + 1;
    jm1 = mod(j-2, N) + 1;
    jm2 = mod(j-3, N) + 1;

    % �Ľ����Ĳ�ֹ�ʽ
    d4u(j) = (u(jm2) - 4*u(jm1) + 6*u(j) - 4*u(jp1) + u(jp2)) / dx^4;
end

end


