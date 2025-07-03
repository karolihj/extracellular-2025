function [S_pattern] = set_up_sparsity_pattern(N, Ns, V_idx)
%[S_pattern] = set_up_sparsity_pattern(N, Ns, V_idx)

S_pattern = kron(speye(N), ones(Ns));
for i=1:N-1
    S_pattern((i-1)*Ns + V_idx,i*Ns + V_idx) = 1;
    S_pattern(i*Ns + V_idx,(i-1)*Ns + V_idx) = 1;
end

% Ensure sparse
S_pattern = sparse(S_pattern);

end