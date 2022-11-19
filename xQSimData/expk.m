function C = Expk(~)
% C = EXPK(); Returns experimental data for the probability
% of clicks per channel. 

load('expk.mat','Nc','Count');
C = Nc;
end 
