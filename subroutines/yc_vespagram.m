function [ves] = yc_vespagram(din,h,s,t,order)
%%% Calculate vespagram
% Written by Yangkang Chen
% Dec, 2018
% 
% order: Order of Nth-root stack (if 1, reverts to linear Radon transform)

nt=length(t);
dt=t(2)-t(1);

Param.h=h;
Param.v=s;
Param.nt=nt;
Param.dt=dt;
Param.type=4;
Param.order=1/order;

ves=yc_radon_pseudo(din,Param,-1);
ves=sign(ves).*abs(ves).^order;

end

