function [d d_eq d_op] = compute_nn_pw_distances(x,y,c)
%
% CALLED BY ANALYZE_SINGLE_MAP.M
%
% USAGE : 
% %FUNCTION  [d d_eq d_op] = compute_nn_pw_distances(x,y,c)
%
%
% INPUT PARAMETERS: 
% x ... x coordinates of all pinwheels
% y ... y coordinates of all pinwheels
% c ... chirality (charge) of all the pinwheels, either 1 for positive or
%       -1 for negative charge
%
%
%
% OUTPUT PARAMETERS:
% d ... nn distances irrespect of which chirality (topological charge)
% d_eq ... nn distances for pinwheels of equal chirality (topological charge)
% d_op ... nn distances for pinwheels of opposite chirality (topological charge)
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	l = length(x);
    for t = 1:l		
            xa = x(t);
            ya = y(t);						
            d(t) = find_nn(x,y,xa,ya);						
    end

    % d same chirality

	ind_p = find(c== 1);
	ind_n = find(c==-1);

	xp = x(ind_p);
	yp = y(ind_p);

	xn = x(ind_n);
	yn = y(ind_n);

	l = length(xp);

    for t = 1:l

            xa = xp(t);
            ya = yp(t);

            dp(t) = find_nn(xp,yp,xa,ya);

    end

	l = length(xn);

    for t = 1:l

            xa = xn(t);
            ya = yn(t);

            dn(t) = find_nn(xn,yn,xa,ya);

    end

	d_eq = [dp dn];
	
    % d opposite chirality

	ind_p = find(c== 1);
	ind_n = find(c==-1);

	xp = x(ind_p);
	yp = y(ind_p);

	xn = x(ind_n);
	yn = y(ind_n);

	l = length(xp);

    for t = 1:l

            xa = xp(t);
            ya = yp(t);

            dp(t) = find_nn(xn,yn,xa,ya);

    end

	l = length(xn);

    for t = 1:l

            xa = xn(t);
            ya = yn(t);

            dn(t) = find_nn(xp,yp,xa,ya);

    end

	d_op = [dp dn];

end %%%% END OF COMPUTE_NN_DISTANCES.M
    


%%%%%%%%%%%%%%%%%%%%%%%%%% LIST OF USED FUNCTIONS %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
function d = find_nn(x,y,xa,ya)
	
	xt = x - xa;
	yt = y - ya;

	r = sqrt(xt.^2 + yt.^2);
	ind = find(r >0);
	d = min(r(ind));
	
end	