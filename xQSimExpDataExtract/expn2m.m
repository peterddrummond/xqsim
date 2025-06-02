function C = expn2m(~)
% C = EXPN2M(); Returns comparison experimental data for all second-order
% correlation moments.

load('expn2m.mat','n2p','Count');
C = n2p;
end 