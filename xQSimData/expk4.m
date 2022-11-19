function C = k4_exp(~)
% C = K4_EXP(); Returns comparison experimental data for the probability
% of total number of counts with four-fold partition. 

load('exp_cp4.mat','Cp4','Count');
C = Cp4;
end 

