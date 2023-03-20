function [an,th0,ph0,sig]=pw_parameter(f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USAGE:
%
% function [an,th0,ph0,sig]=pw_parameter(f)
% This function computes chirality, overrepresented orientations and
% anisotropy parameter of a pinwheel
%
% INPUT PARAMETER:
% f ... 4x1 vector containing the local directional derivateds ([grad_x Re(z),grad_y Re(z), grad_x Im(z), grad_y Im(z)])
%
% OUTPUT PARAMETERS: 
% an  ... pinwheel anisotropy 
% th ...  overrepresented orientation givein i(-Pi/2,Pi/2)
% ph0 ...    overrepresented orientation give in (0,2 Pi)
% sig ... chirality +/- 1
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a=f(1,:); b=f(2,:); c=f(3,:); d=f(4,:);
    M   = [a,b;c,d];

    al  = c^2+d^2;
    ga  = a^2+b^2;
    be  = -(a*c+b*d);

    an  = 1-sqrt( (al-sqrt(4*be^2+(al-ga)^2)+ga) / (al+sqrt(4*be^2+(al-ga)^2)+ga));
    th0 = -atan( 2*be/(ga-al+sqrt(4*be^2+(al-ga)^2)) );
    vec = M\[2*be^2;be*(ga-al-sqrt(4*be^2+(al-ga)^2))];%involves inversion of matrix M!
    com = vec(1)+1i*vec(2);
    ph0 = angle(com);
    sig = sign(det(M));

% 
%     rho=atan(c/a);
%     sigma=-atan(b/d)-rho;
%     amp=a/cos(rho);
%     alp=-c/(amp*sin(rho+sigma));

end