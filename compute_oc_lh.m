function OC_LH = compute_oc_lh(FW_LH,BK_LH,length_Oin)
% function to compute occupation likelihoods
OC_LH = zeros(3,length_Oin);

for t=1:length_Oin
    for i=1:3
        OC_LH_denominator = 0;
        for j=1:3
            OC_LH_denominator = OC_LH_denominator + (FW_LH(j,t)*BK_LH(j,t));
        end
        OC_LH(i,t) = (FW_LH(i,t)*BK_LH(i,t))/OC_LH_denominator;
    end
end
return;