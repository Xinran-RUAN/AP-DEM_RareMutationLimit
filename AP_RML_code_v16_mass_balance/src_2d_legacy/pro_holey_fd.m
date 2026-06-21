function d = pro_holey_fd(p,cc1,cc2,zeta)

% parameters:
%     input, real P, one or more points.
%     output, real D, the signed distance of each point to the boundary of
%     the region.

% g1 = dcircle (p, 0, 0, 1);
% g2 = drectangle ( protate ( p, pi/12), -1 , 1, 0, 1);
% g3 = drectangle ( protate ( p ,-pi/12), -1, 1, -1, 0);
% g4 = drectangle ( protate ( pshift ( p, -0.9, 0), -pi/4), 0, 0.2, 0, 0.2);
% g5 = dcircle ( p, 0.6, 0, 0.1);
% d = ddiff( ddiff( ddiff( ddiff (g1, g2), g3), g4), g5);

g1 = dcircle (p, 0, 0, 1);
g2 = dcircle (p, cc1, 0, zeta);
g3 = dcircle (p, cc2, 0, zeta);

d = ddiff( ddiff (g1, g2), g3);

return

end

