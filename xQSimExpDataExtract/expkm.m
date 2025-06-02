function C = km_exp(~)
% C = KM_EXP(); Returns comparison experimental data for the probability
% of up to K-th moment. 

load('exp_km.mat','Kp','Count');
C = Kp;
end 
