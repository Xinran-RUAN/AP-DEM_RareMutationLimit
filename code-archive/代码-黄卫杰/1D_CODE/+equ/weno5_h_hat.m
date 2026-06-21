function grad_sq = weno5_h_hat(u, dx)

pL = equ.WENO5_left_RXR(u, dx);
pR = equ.WENO5_right_RXR(u, dx);
grad_sq = equ.RF_flux(pR, pL);

end

