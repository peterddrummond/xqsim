function C = expk3ms(~)

% C = EXPK3MS(); Returns comparison of experimental data a subset of all
% possible combiantions of output modes for the probability of 3nd-order
% moments. Product is taken over sequential modes.

load('expk3ms.mat','K3sp','Count');
C = K3sp;
end 