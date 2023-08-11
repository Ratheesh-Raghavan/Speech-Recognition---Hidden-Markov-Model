%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                    %
%                                                                       %
% MATLAB Script.                                                        %
% Function: Create a Hidden Marcov Model and perform Forward-Backward   %
%           procedure and use Baum-Welch equations to re-estimate the   %
%           model parameters (train) based on a a set of observation    %
%           sequences.                                                  %
%           Number of possible states = 3 (1,2,3)                       %
%           Number of possibe observations = 8 (1,2,3,4,5,6,7,8)        %
%                                                                       %
% Written by: Ratheesh Raghavan                                         %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

% Hidden Marcov Model
A = [0.80 0.05 0.05; 0.15 0.70 0.05; 0.05 0.10 0.75];
B = [0.30 0.20 0.25 0.10 0.05 0.05 0.00 0.05;...
     0.00 0.00 0.05 0.05 0.05 0.20 0.30 0.35;...
     0.00 0.05 0.05 0.20 0.35 0.25 0.05 0.05];
pi = [0.40 0.25 0.35];
eta = [0.1; 0.1; 0.1];

% Observation Sequences 
O{1} = [6 4 2 1 3];
O{2} = [8 7 6 8 7 3 6 3];
O{3} = [5 4 6 7];
O{4} = [5 7 5 5 8 8 6]; 
O{5} = [3 2 2 2 2 1 8 7 7 3];

% Containers to hold forward, backword, occ and trans likelihoods
FW_LH = [];
BK_LH = [];
OC_LH = [];
TR_LH = {};

% Initialize the accumulators
accum_estB_numerator = zeros(3,8);
accum_OC_LH = zeros(3,1);
accum_TR_LH = zeros(3,3);
accum_OC_LH_t_1 = zeros(3,1);
accum_pi_numerator = zeros(1,3);
accum_eta_numerator = zeros(3,1);

% Looping for each observation sequence 
for i=1:length(O)
% Compute Forward Likelihoods and Overall Likelihood
    FW_LH = compute_fw_lh(A,B,pi,O{i});
    fprintf('O%d - Forward Likelihoods\n',i);
    disp(FW_LH); 
    Overall_LikelihoodFW = FW_LH(1,length(O{i}))*eta(1,1) +...
                           FW_LH(2,length(O{i}))*eta(2,1) +...
                           FW_LH(3,length(O{i}))*eta(3,1);
    fprintf('O%d - Overall Likelihood (forward procedure) = %f\n\n',...
            i,Overall_LikelihoodFW);           
% Compute Backward Likelihoods
    BK_LH = compute_bk_lh(A,B,eta,O{i});
    fprintf('O%d - Backward Likelihoods\n',i);
    disp(BK_LH); 
    Overall_LikelihoodBK = BK_LH(1,1)*pi(1,1)*B(1,(O{i}(1))) +...
                           BK_LH(2,1)*pi(1,2)*B(2,(O{i}(1))) +...
                           BK_LH(3,1)*pi(1,3)*B(3,(O{i}(1)));
    fprintf('O%d - Overall Likelihood (backward procedure) = %f\n\n',...
            i,Overall_LikelihoodBK);
% Compute Occupation Likelihoods
    OC_LH = compute_oc_lh(FW_LH,BK_LH,length(O{i}));
    fprintf('O%d - Occupation Likelihoods\n',i);
    disp(OC_LH);        
% Compute Transition Likelihoods
    TR_LH = compute_tr_lh(FW_LH,BK_LH,A,B,O{i});
    for tr=1:length(O{i})-1
        fprintf('O%d - Transition Likelihoods ',i);
        fprintf('for transition from t=%d to t=%d\n',tr,tr+1);
        disp(TR_LH{tr}); 
    end
% Accumulate the Occupation and Transition Likelihoods
    for t=1:length(O{i})
        for s=1:3
            accum_estB_numerator(s,O{i}(t)) = ... 
            accum_estB_numerator(s,O{i}(t)) + OC_LH(s,t);
        end
    end
    for s=1:3
        accum_OC_LH(s,1) = accum_OC_LH(s,1) + sum(OC_LH(s,1:length(O{i})));
        accum_OC_LH_t_1(s,1) = accum_OC_LH_t_1(s,1) + ...
                               sum(OC_LH(s,1:length(O{i})-1)); 
        accum_pi_numerator(1,s) = accum_pi_numerator(1,s) + OC_LH(s,1);
        accum_eta_numerator(s,1) = accum_eta_numerator(s,1) +...
                                   OC_LH(s,length(O{i}));
    end
    for t=1:length(O{i})-1
        accum_TR_LH = accum_TR_LH + TR_LH{t};
    end
end

% Re-estimate State-Transition Probabilities
A_bw = accum_TR_LH ./accum_OC_LH_t_1;
disp('Estimated new state transition probabilities A_bw');
disp(A_bw);
pi_bw = accum_pi_numerator ./length(O);
disp('Estimated new state transition probabilities pi_bw');
disp(pi_bw);
eta_bw = accum_eta_numerator ./length(O);
disp('Estimated new state transition probabilities eta_bw');
disp(transpose(eta_bw));

% Re-estimate Output Probabilities
B_bw = accum_estB_numerator ./accum_OC_LH;
disp('Estimated new output probabilities B_bw');
disp(B_bw);

% Plot initial and re-estimated Output Probabilities
figure(1);
X = 1:1:8;
subplot(1,3,1);
Y = B(1,:);
plot(X,Y,'r');
title("State 1");
xlabel("Observation Type (1 - 8)");
ylabel("Output Probability");
hold on;
subplot(1,3,2);
Y=B(2,:);
plot(X,Y,'r');
title("State 2");
xlabel("Observation Type (1 - 8)");
ylabel("Output Probability");
hold on;
subplot(1,3,3);
Y=B(3,:);
plot(X,Y,'r');
title("State 3");
xlabel("Observation Type (1 - 8)");
ylabel("Output Probability");
hold on;
subplot(1,3,1);
Y = B_bw(1,:);
plot(X,Y,'b');
subplot(1,3,2);
Y=B_bw(2,:);
plot(X,Y,'b');
subplot(1,3,3);
Y=B_bw(3,:);
plot(X,Y,'b');
sgtitle("Initial and Re-estimated Output Probabilities - Comparison");
hold off;

% Compuet most likely path for each observation sequence
% (this is for comparison with Viterby reference training results)
fprintf('Most Likely Paths based on Re-estimated HMM parameters\n');
for i=1:length(O)
    FW_LH = compute_fw_lh(A_bw,B_bw,pi_bw,O{i});
    BK_LH = compute_bk_lh(A_bw,B_bw,eta_bw,O{i});
    OC_LH = compute_oc_lh(FW_LH,BK_LH,length(O{i}));
    % Compute most likely path
    [MaxVal,MaxState] = max(OC_LH(:,1:length(O{i})));
    disp(MaxState);
end
