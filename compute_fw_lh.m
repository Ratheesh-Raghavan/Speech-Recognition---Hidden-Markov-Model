function FW_LH = compute_fw_lh(A,B,pi,Oin)
% function to compute forward likelihoods
FW_LH = zeros(3,length(Oin));
for i=1:3
    FW_LH(i,1) = pi(1,i)*B(i,Oin(1)); 
end

for t=2:length(Oin)
    for i=1:3
        tot_tran_in_prob = 0;
        for j=1:3
            tran_in_prob = FW_LH(j,t-1)*A(j,i);
            tot_tran_in_prob = tot_tran_in_prob + tran_in_prob;
        end
        FW_LH(i,t) = B(i,Oin(t)) * tot_tran_in_prob;
    end
end
return;