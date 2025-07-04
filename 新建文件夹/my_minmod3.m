% u is row vector 
function[mmu] = my_minmod3(u1, u2, u3)

mmu = zeros(size(u1));
idx_pos = (u1 > 0) & (u2 > 0) & (u3 > 0);
idx_neg = (u1 < 0) & (u2 < 0) & (u3 < 0);
Up = cat(1, u1(idx_pos), u2(idx_pos), u3(idx_pos));
mmu(idx_pos) = min(Up, [], 1);
Un = cat(1, u1(idx_neg), u2(idx_neg), u3(idx_neg));
mmu(idx_neg) = max(Un, [], 1);
% mm is zero elsewhere (where signs differ or either is zero)
