function C = expk2ms(~)

% C = EXPK2MS(); Returns comparison of experimental data a subset of all
% possible combinations of output modes for the probability of 2nd-order
% moments. Product is taken over sequential modes e.g. data point one is
% click product for modes 1,2, data point two is click product for modes
% 2,3 etc.


load('expk2ms.mat','K2sp','Count');
C = K2sp;
end 