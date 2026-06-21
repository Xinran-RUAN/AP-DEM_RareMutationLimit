x = 0:0.01:1-0.01;
y = sin(2*pi*x);
dy_l = WENO5_left(y,0.01);
dy_r = WENO5_right(y,0.01);
dy = 2 * pi * cos(2*pi*x);
figure(1)
plot(x, dy_l - dy); hold on;
plot(x, dy_r - dy, 'r--'); hold off; % dy_l 与 dy_r分别为单侧信息，加和之后才是对应dy
