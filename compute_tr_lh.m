function TR_LH = compute_tr_lh(FW_LH,BK_LH,A,B,Oin)
% function to compute transition likelihoods
TR_LH = {};

for t=1:length(Oin)-1
    TR_LH_denom = 0;
    for k=1:3
        for w=1:3
            TR_LH_denom = TR_LH_denom +...
                          FW_LH(k,t)*A(k,w)*BK_LH(w,t+1)*B(w,Oin(t+1));
        end
    end
    for i=1:3
        for j=1:3
            TR_LH{t}(i,j) = FW_LH(i,t)*A(i,j)*BK_LH(j,t+1)*B(j,Oin(t+1))...
                            /TR_LH_denom;
        end
    end
end
return;