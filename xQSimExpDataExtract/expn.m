function C = expn(~)
% C = EXPN(); Returns experimental data for the probability
% of photons per channel. 

load('expn.mat','Nc','Count');
C = Nc;
end 
