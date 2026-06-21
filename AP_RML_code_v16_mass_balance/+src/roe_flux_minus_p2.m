function Hhat = roe_flux_minus_p2(pPlus, pMinus)
%ROE_FLUX_MINUS_P2 Roe/entropy-fix flux for H(p) = -p^2.
H = @(p) -p.^2;
Hp = @(p) -2*p;
pPlus = pPlus(:).';
pMinus = pMinus(:).';
pbar = 0.5*(pPlus + pMinus);
pmin = min(pMinus, pPlus);
pmax = max(pMinus, pPlus);
noChange = (Hp(pmin) >= 0 & Hp(pmax) >= 0) | (Hp(pmin) <= 0 & Hp(pmax) <= 0);
pStar = pPlus;
pStar(Hp(pbar) >= 0) = pMinus(Hp(pbar) >= 0);
Hhat = zeros(size(pPlus));
Hhat(noChange) = H(pStar(noChange));
alpha = 2*max(abs(pPlus), abs(pMinus));
Hhat(~noChange) = H(0.5*(pPlus(~noChange) + pMinus(~noChange))) ...
                  - 0.5*alpha(~noChange).*(pPlus(~noChange)-pMinus(~noChange));
end
