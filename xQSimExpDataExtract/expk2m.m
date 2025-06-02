function C = expk2m(~)
% C = EXPK2M(); Returns comparison experimental data for all second-order
% correlations.

load('expk2m.mat','K2p','Count');
C = K2p;
end 