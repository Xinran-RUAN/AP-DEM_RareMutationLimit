function grad_sq = weno5_h_hat(u, dx)

pL = WENO5_left_RXR(u, dx);
pR = WENO5_right_RXR(u, dx);
grad_sq = RF_flux(pR, pL);

end

