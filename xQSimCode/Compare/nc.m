function C = nc(p)
% C = NC(p); generates comparison numbers of  clicks per channel
% assumes output unchanged from input, apart from transmission factor t
% eg, identity matrix or uniform thermal+unitary

C = (sinh(p.r').*p.t').^2;
end                                         %End nc function