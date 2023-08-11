function BK_LH = compute_bk_lh(A,B,eta,Oin)
% function to compute backward likelihoods
BK_LH = zeros(3,length(Oin));
for i=1:3
    BK_LH(i,length(Oin)) = eta(i,1); 
end

for t=length(Oin)-1:-1:1
    for i=1:3
        for j=1:3
            BK_LH(i,t) = BK_LH(i,t) + BK_LH(j,t+1)*A(i,j)*B(j,(Oin(t+1)));
        end
    end
end
return;