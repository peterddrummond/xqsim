function C = expk3m(~)
% C = EXPK3M(); Returns comparison experimental data for all third-order
% correlations.

load('expk3m.mat','K3p','Count');
C = K3p;
end 