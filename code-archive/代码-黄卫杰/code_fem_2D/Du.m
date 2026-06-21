function[du] = Du(u)
theta = 1;
u_l = [u(end), u(1:end-1)];
u_r = [u(2:end), u(1)];
du_l = u - u_l;
du_r = u_r - u;
du = my_minmod3(theta * du_l, theta * du_r, 0.5 * (du_l + du_r));
