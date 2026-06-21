% 显式部分
function rhs = RHS(u, dx, H)
    flux = monotone_flux(u, dx);
    rhs = flux - H;
end
