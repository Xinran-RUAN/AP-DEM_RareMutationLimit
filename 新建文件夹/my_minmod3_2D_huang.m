% W is a matrix 
function[mmW] = my_minmod3_2D_huang(W1, W2, W3)
[m, n] = size(W1);
W1_vector = reshape(W1, [], 1);
W2_vector = reshape(W2, [], 1);
W3_vector = reshape(W3, [], 1);
mmW_vector = zeros(size(W1_vector));
idx_pos = (W1_vector > 0) & (W2_vector > 0) & (W3_vector > 0);
idx_neg = (W1_vector < 0) & (W2_vector < 0) & (W3_vector < 0);
Up = cat(2, W1_vector(idx_pos), W2_vector(idx_pos), W3_vector(idx_pos));
mmW_vector(idx_pos) = min(Up, [], 2);
Un = cat(2, W1_vector(idx_neg), W2_vector(idx_neg), W3_vector(idx_neg));
mmW_vector(idx_neg) = max(Un, [], 2); 
% mm is zero elsewhere (where signs differ or either is zero)
mmW = reshape(mmW_vector, m, n);    
