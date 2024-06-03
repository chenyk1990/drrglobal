function [out] = yc_radon_pseudo(in,Param,operator);
%Forward and Adjoint operators for Radon RT (Pseudo) in the time domain.
%  IN   	in:   		intput data 
%       	Param:  		parameters combination
%		Param.h:		offset
%		Param.v:		velocity
%		Param.nt:		number of samples
%		Param.dt:		time interval
%		Param.type:    1: linear 2: parablic 3: hyperbolic 4:v->slowness
%
%	   	operator: 
%% 			operator =  1 means impute is m(tau,v) and output is d(t,x) FORWARD  OP
%% 			operator = -1 means input is d(t,x) and output is m(tau,v)  ADJOINT  OP 
%      
%  OUT   out:  		output data
% 
%  Copyright (C) 2018 Zhejiang University
%  Copyright (C) 2018 Yangkang Chen
%
% Dot test example (in the case vespagram) not succeed:
%  nt=500;
%  nv=100;
%  nh=100;
%  h=1:nh;
%  dt=1;
%  type=1;
%  
%  Param.h=h;
%  Param.nt=nt;
%  Param.dt=dt;
%  Param.v=linspace(-5,10,nv);
%  Param.type=type;
%  Param.order=1;
% 
%  m1 = rand(nt,nv); 
% [d1 ] = radon_coda(m1,Param,1);
% 
%  d2 = rand(nt,nh); 
% [m2 ] = radon_coda(d2,Param,-1);
% 
% dot1 = sum(sum(d1.*d2))
% dot2 = sum(sum(m1.*m2))


%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/


h = Param.h;
v = Param.v;
nt = Param.nt;
dt = Param.dt;
type=Param.type;

if isfield(Param,'order')
   order=Param.order;
else
   order=1;
end


 nh = length(h);
 nv = length(v);

 if operator == -1; m = zeros(nt,nv); end;
 if operator ==  1; d = zeros(nt,nh); end;

 if operator == -1; d = in; end; 
 if operator ==  1; m = in; end; 

 hmax=max(abs(h));
 
  for itau = 1:nt
 
    for ih = 1:nh
 
      for iv = 1:nv

%% This can also be replaced by Parabolic or linear integration 

    switch type
        case 1
             t = (itau-1)*dt + h(ih)/v(iv);	   	
            it = floor(t/dt)+1;            
        case 2
             t = (itau-1)*dt + h(ih)*h(ih)*v(iv)/hmax/hmax;   %curvature
            it = floor(t/dt)+1; 
            %if(it<=0) it=1;end
        case 3
            t = sqrt (((itau-1)*dt)^2 + (h(ih)/v(iv))^2 ) ;
            it = floor(t/dt)+1; 
        case 4
             t = (itau-1)*dt + h(ih)*v(iv);	 %here v means slowness  	
            it = floor(t/dt)+1;     
        otherwise
             t = sqrt (((itau-1)*dt)^2 + (h(ih)/v(iv))^2 ) ;
            it = floor(t/dt)+1;   
    end

%% This code does not apply interpolation 

 if (it<=nt && it>0);
	if operator==-1;     m(itau,iv) = m(itau,iv) + sign(d(it,  ih))*abs(d(it,  ih)).^(order)/nh;  end
	if operator== 1;     d(it,  ih)   = d(it,ih) + sign(m(itau,iv))*abs(m(itau,iv)).^(1/order)/nh; end
% 	if operator==-1;     m(itau,iv) = m(itau,iv) + d(it,  ih).^(order)/nh;  end
% 	if operator== 1;     d(it,  ih)   = d(it,ih) + m(itau,iv).^(order)/nh; end
 end

  end
  end
  end

 if operator == 1; out = d; end; 
 if operator ==-1; out = m; end; 

 

