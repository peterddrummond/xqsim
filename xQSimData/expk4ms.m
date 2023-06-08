function C = expk4ms(~)

% C = EXPK4MS(); Returns comparison of experimental data a subset of all
% possible combiantions of output modes for the probability of 4th-order
% moments. Product is taken over sequential modes.

load('expk4ms.mat','K4sp','Count');
C = K4sp;
end 