function[dW] = DW(W)
theta = 1;
W_l = [W(:, end), W(:, 1:end-1)];
W_r = [W(:, 2:end), W(:, 1)];
dW_l = W - W_l;   
dW_r = W_r - W;
dW = my_minmod3_2D_huang(theta * dW_l, theta * dW_r, 0.5 * (dW_l + dW_r));