function P = introduce_variation(P, G)
%P = introduce_variation(P, G)

idx = [1,2,6,7,8,9,10,11,12,18,22,23,25]; % Parameter indices

load(sprintf('random_numbers/N%d.mat', G.Nm), 'a')

P(idx,:) = P(idx,:).*(1 + 0.5*G.vary_percentage*(-1+2*a));

end