function [ D1 ] = fk_pocs(D,D0,mask,thr,niter,type)
%FK POCS: FK POCS interpolation (POCS + IST)
%  IN   D:   	intput data 
%       D0:   	initial data
%               mask: masking operator
%               thr: percentile threshold
%               niter: number of iterations
%               type: POCS or IST
%               'pocs'-POCS and 'ist'-IST
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
D1=D0;
 switch type 
     case 'pocs'
         for iter=1:niter
            D1=D.*mask+(1-mask).*yc_fkt(D1,'ps',thr); 
         end
         
     case 'ist'
         for iter=1:niter
           Dtmp=D1+D-mask.*D1;  
           D1=yc_fkt(Dtmp,'ps',thr);  
         end
         
  otherwise 
    error('Invalid argument value.')      
    
 end
return
