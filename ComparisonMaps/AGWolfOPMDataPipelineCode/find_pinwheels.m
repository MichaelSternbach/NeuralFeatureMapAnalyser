function [count,aniso,x_angle,PWxList,PWyList,signList, contours] = find_pinwheels(z,rm,roi_matrix,verboseMode)
%
% CALLED BY ESTIMATE_LOCAL_PW_DENSITY.M
%
% USAGE:
%    function [count,aniso,x_angle,PWxList,PWyList,signList] = pw_finder(z,rm,roi_matrix,verboseMode)
%
% DESCRIPTION:
%
% This function finds all pinwheels in a complex-valued matrix z by
% computing intersections between zero-line contours of real and imaginary
% part.
% 
%
% INPUT PARAMETERS:
% z         ... complex-valued 2D matrix
% rm        ... radius for pinwheel locations to be exclude around the
%               edges of the ROI (rm = 0 in Kaschube et al. 2010, Keil et al. 2011, Science)
% roi_matrix... 2D matrix, containing ones inside the ROI, and zeros otherwise
%
% verboseMOde... flag,if 0 faster processing is achieved, if 1 maps and contourlines
%                are displayed for cross-validation
%
%
% OUTPUT PARAMETERS:
% count     ... number of pinwheels in the map
% aniso     ... anisotropy parameter for each pinwheel
% x_angle   ...
% PWxList   ...
% PWyList   ...
% signList  ...
% contours  ... structure with the contours and their indices
%
%
%
%  Copyright (C) 2014, Max-Planck-Institut for Dynamics and Self-organization, The  
%  Nonlinear Dynamics Group. This software may be used, copied, or 
%  redistributed as long as it is not sold and this copyright notice is 
%  reproduced on each copy made. This routine is provided as is without 
%  any express or implied warranties whatsoever.  
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    warning off MATLAB:divideByZero; %#ok<*RMWRN>
    count=0;
    x_angle=[];
    aniso =[];
    nx =2;
    ny =2;
    
    % input checks
    if nargin <2
        rm      = 0 ;  %radius in pixel, distance from roi_border until which pinwheels are 
                       %discarded, rm=0 is a good choice usually if data doesn't taper off at edges
    end

    [sy,sx]  =size(z);

    if nargin < 3
        roi_matrix = ones(size(z));
    end
    if nargin <4
        verboseMode = false;
    end

    % initialize the pinwheel position and chirality lists
    PWxList=[];
    PWyList=[];
    signList=[];

    % These will be variables to avoid double counting of pinwheels
    doppeltestx=0;   
    doppeltesty=0;
    doppeldet=0;

    %Compute gradient field of z, will be used to determine chirality of
    %pinwheels
    [grx,gry]=gradient(real(z));
    [gix,giy]=gradient(imag(z));

    % Visualize contours, if desired, this will slow down pw finding
    % considerably, so set verboseMode to 0 in application that require
    % speed
    if(verboseMode)
        figure();
        subplot(nx,ny,1);
        contour(real(z),[0 0]);
        axis image,axis off;
        axis xy;
        drawnow;


        co=hsv(1024);
        co(1024,:)=[1 1 1];
        colormap(co)

        subplot(nx,ny,2)
        contour(imag(z),[0 0]);
        axis image,axis off;
        axis xy;
        drawnow;

        subplot(nx,ny,4),hold on
        contour(real(z),[0 0],'k');
        contour(imag(z),[0 0],'k');
        axis image,axis off, axis xy,hold on;


        subplot(nx,ny,3)
        imagesc(angle(z));
        axis image,axis off, axis xy,hold on
        contour(real(z),[0 0],'k');
        contour(imag(z),[0 0],'w');
        drawnow


        fig2 = figure();
        imagesc(angle(z));
        axis image,axis off, axis xy,hold on
        contour(real(z),[0 0],'k');
        contour(imag(z),[0 0],'w');
        drawnow
        colormap(hsv);
    end

    
    % Calculate Contours - Both Real and Imaginary components
    [cReal cRealIndex cRealElements ] = estimateContour(real(z));
    [cImag cImagIndex cImagElements ] = estimateContour(imag(z));

    
    % Eliminate contours which are less than three points long
    cRealIndex = cRealIndex(cRealElements > 2);
    lre = sum(cRealElements > 2);
    cRealElements = cRealElements(cRealElements > 2);

    cImagIndex = cImagIndex(cImagElements > 2);
    lim = sum(cImagElements > 2);
    cImagElements = cImagElements(cImagElements > 2);


    % calculate distances btw. contours & rectangles which contain the contours in advance!
    % 1. contours of imaginary part
     for s= 1:lim

         currentList = (cImagIndex(s)+1):(cImagIndex(s)+cImagElements(s));
         xim = cImag(1,currentList);
         yim = cImag(2,currentList);
         len=cImagElements(s); %no elements in the contour

         Xim(s,1:len) =xim(1:len);
         Yim(s,1:len) =yim(1:len);  % cut off the NaNs
         Nim(s) = len;

         dist2=max((xim(1:len-1)-xim(2:len)).^2+(yim(1:len-1)-yim(2:len)).^2);%% largest spacing between two contour points
         Dist2(s)=dist2;

         im_4(s,1:4)= [min(xim) max(xim) min(yim) max(yim)];

         %%% uncomment, for visualizing the rectangles
         %figure(3)
         %rectangle('position',[min(xim),min(yim),max(xim)-min(xim),max(yim)-min(yim)],'edgecolor','r')
         %hold on, drawnow
         

     end % for s= 1:lim


     % 2. contours of imaginary part
     for s= 1:lre

         currentList = (cRealIndex(s)+1):(cRealIndex(s)+cRealElements(s));
         xre = cReal(1,currentList);
         yre = cReal(2,currentList);
         len=cRealElements(s);

         Xre(s,1:len) =xre(1:len);
         Yre(s,1:len) =yre(1:len);  % cut off the NaNs
         Nre(s) = len; % number of elements in the contour

         dist1=max((xre(1:len-1)-xre(2:len)).^2+(yre(1:len-1)-yre(2:len)).^2); %% largest spacing between two contour points
         Dist1(s)=dist1;

         re_4(s,1:4)= [min(xre) max(xre) min(yre) max(yre)];

         %%% uncomment, for visualizing the rectangles
         %figure(3)
         %rectangle('position',[min(xre),min(yre),max(xre)-min(xre),max(yre)-min(yre)],'edgecolor','k')
         %hold on, drawnow

     end % for s= 1:lre


     % now find out which contours meet, i.e. which of the rectangles overlap
     M = zeros(lre,lim);

     for t =1:lre

         min_x_re = re_4(t,1);
         max_x_re = re_4(t,2);
         min_y_re = re_4(t,3);
         max_y_re = re_4(t,4);


         for s=1:lim

             min_x_im = im_4(s,1);
             max_x_im = im_4(s,2);
             min_y_im = im_4(s,3);
             max_y_im = im_4(s,4);


             % if this expression is true, the contours meet! 
             if ( (((min_x_re <= min_x_im) && (min_x_im <= max_x_re)) || ...
                   ((min_x_re <= max_x_im) && (max_x_im <= max_x_re)) || ...
                   ((min_x_im <= min_x_re) && (min_x_re <= max_x_im)) || ...
                   ((min_x_im <= max_x_re) && (max_x_re <= max_x_im))) && ...
                  (((min_y_re <= min_y_im) && (min_y_im <= max_y_re)) || ...
                   ((min_y_re <= max_y_im) && (max_y_im <= max_y_re)) || ...
                   ((min_y_im <= min_y_re) && (min_y_re <= max_y_im)) || ...
                   ((min_y_im <= max_y_re) && (max_y_re <= max_y_im))) )

                 % set the 'meet' matrix M to one for this entry
                 M(t,s)=1;
             end
         end
     end

     
    if(verboseMode) 
        figure(1)
    end
    
    % Now go through all the contours, and ask for the contours they meet
    % whether they actually cross
    
    % Loop over all contours of the real part
    for t=1:lre

        lr=Nre(t);
        xre = Xre(t,:);
        yre = Yre(t,:);
        xre = xre(1:lr);
        yre = yre(1:lr);

        dist1=Dist1(t);

        mm  = M(t,:);
        s2t = find(mm==1);
        
        % Loop over all contours of the imaginary part that contour no. t meets
        for ss= 1:length(s2t);
            
            s = s2t(ss);
            

            li   = Nim(s);

            xim = Xim(s,:);
            yim = Yim(s,:);
            xim = xim(1:li);
            yim = yim(1:li);

            dist2=Dist2(s);

            de = 2*max(dist1,dist2);
            
            
            % loop over all points of the real part contour
            for p=1:lr-1

                % find points on imaginary part contour (that meets the real part contour) that come close to the point xre(p) of the real part contour
                dist = (xre(p)-xim).^2 + (yre(p)-yim).^2;

                index= find(dist<=de); % find points of the imaginary part contour that come close to the point xre(p) of the real contour

                % go through all of those points and look for the crossing 
                for q=1:length(index);

                    if index(q) < li;

                        cx=xim(index(q));
                        cy=yim(index(q));
                        %%figure(2),plot(cx,cy,'dm'),figure(1)


                        ax=xre(p);
                        ay=yre(p);
                        bx=xre(p+1);
                        by=yre(p+1);
                        dx=xim(index(q)+1);
                        dy=yim(index(q)+1);
                        
                        [answer,m,n] =  dotheycross(ax,ay,bx,by,cx,cy,dx,dy);

                        if answer
                            % Crossin found!
                            %%%% check whether there is direct backcrossing
                            %%%% in the next contour point, if so, ignor
                            %%%% the crossing!
                            answer2 = 0;
                            answer3 = 0;
                            answer4 = 0;
                            answer5 = 0;
                            answer6 = 0;
                            
                            if p < lr-1
                                a2x = xre(p+1);
                                a2y = yre(p+1);
                                b2x = xre(p+2);
                                b2y = yre(p+2);
                                
                                answer2 =  dotheycross(a2x,a2y,b2x,b2y,cx,cy,dx,dy);
                                
                                if index(q) < li-1
                                    c4x = xim(index(q)+1);
                                    c4y = yim(index(q)+1);
                                    d4x = xim(index(q)+2);
                                    d4y = yim(index(q)+2);
                                    answer5 =  dotheycross(a2x,a2y,b2x,b2y,c4x,c4y,d4x,d4y);
                                    %%%% Both contours have one point on
                                    %%%% other side in this case
                                end
                            end
                            
                            %%% Go one point back, if possible and check if
                            %%% crossing!
                            if p > 1
                                a3x = xre(p-1);
                                a3y = yre(p-1);
                                b3x = xre(p);
                                b3y = yre(p);
                                
                                answer3 =  dotheycross(a3x,a3y,b3x,b3y,cx,cy,dx,dy);
                                
                                if index(q) > 1
                                    
                                    c6x = xim(index(q)-1);
                                    c6y = yim(index(q)-1);
                                    d6x = xim(index(q));
                                    d6y = yim(index(q));
                                   
                                    answer6 =  dotheycross(a3x,a3y,b3x,b3y,c6x,c6y,d6x,d6y);
                                    %%%% if answer6==1, both contours have one point on
                                    %%%% other side
                                end
                                
                            end
                            %%% Backcrossing of imaginary part 
                            if index(q) < li-1
                                c4x = xim(index(q)+1);
                                c4y = yim(index(q)+1);
                                d4x = xim(index(q)+2);
                                d4y = yim(index(q)+2);
                                
                                answer4 =  dotheycross(ax,ay,bx,by,c4x,c4y,d4x,d4y);
                            end
                            
                            
                            if (answer2==0) && (answer3==0) && (answer4==0) && (answer5==0) && (answer6==0)
                            % all possibilities of direct backcrossing excluded
                            % actual crossing, i.e. pinwheel found
                            % Now determine chirality of pinwheel using
                            % determinant of local gradient field
                                
                                count=count+1;

                                % Compute the actual crossing point,i.e. pinwheel position by
                                % linearly interpolating the contour
                                pinx=cx+m*(dx-cx);
                                piny=cy+m*(dy-cy);
                                PWxList=[PWxList pinx]; % add to list
                                PWyList=[PWyList piny];
                                
                                x=round(pinx); % rounded position value of the pinwheel
                                y=round(piny); % used later to check whether we already have that pinwheel in the list!

                                % transpose coordinates (x->y, y-->x) because of MATLAB
                                % weirdness in output of coordinates 
                                qx=x;
                                qy=y;

                                x=qy;
                                y=qx;

                                a=grx(x,y); %recall: [grx,gry]=gradient(real(z));
                                b=gry(x,y); %recall: [grx,gry]=gradient(real(z));
                                c=gix(x,y); %recall: [gix,giy]=gradient(imag(z));
                                d=giy(x,y); %recall: [gix,giy]=gradient(imag(z));


                                det=a*d-b*c;  %Determinant determines chirality of the pinwheel
                                signList=[signList det];

                                %back transpose coordinates
                                x=qx;
                                y=qy;

                                %check whether pw surrounding does touch the roi_matrix
                                touch_roi_matrix=0;

                                for  phi = (0:.01:2*pi-.01);

                                    xr=round(x+rm*cos(phi));
                                    yr=round(y+rm*sin(phi));

                                    if xr>1 && xr <sx && yr>1 && yr <sy
                                        if roi_matrix(yr,xr)==0
                                            touch_roi_matrix=1;
                                            break;
                                        end;

                                    else touch_roi_matrix=1; break %% near matrix edge!

                                    end
                                end

                                % Check whether we already have the
                                % pinwheel in the list
                                % pinwheel is defined to already be in the
                                % list if:
                                % 1. its position value, rounded to integer, is
                                % equal to another rounded position value found in the previous loops
                                % 2. its chirality is the same as the one
                                % of that previously found pinwheel
                                if   (qx==doppeltestx) && (qy==doppeltesty) && (sign(det)==doppeldet)


                                    count = count-1;
                                    PWxList=PWxList(1:length(PWxList)-1);
                                    PWyList=PWyList(1:length(PWyList)-1);
                                    signList=signList(1:length(signList)-1);
                                else
                                    if(verboseMode)
                                        subplot(nx,ny,4)
                                        hold on;
                                    end


                                    if x<sx && y <sy && x>0 && y>0
                                        if (roi_matrix(y,x)==1) && (touch_roi_matrix==0)
                                            % pinwheel is within roi and not
                                            % touching any borders
                                            if(verboseMode)
                                                if sign(det) > 1
                                                    plot(x,y,'bo');
                                                else
                                                    plot(x,y,'go');
                                                end
                                                axis xy;
                                            end

                                            % Computes pinwheel anisotropy
                                            %an=pw_parameter([a,b,c,d]');
                                            %aniso=[aniso an];
                                            %xangle=acos((a*c+b*d)/(sqrt(a^2+b^2)*sqrt(c^2+d^2)))/pi*180;
                                            %x_angle=[x_angle xangle];

                                        else
                                            % pinwheel is too close to the ROI border, delete from list                                         
                                            %%%% Plot pinwheels at the borders (eliminated) as red stars    
                                            if(verboseMode)
                                                plot(x,y,'r*'),axis xy;
                                            end

                                            %%%% Elimnate the pinwheels at the borders of the ROI    
                                            count = count-1;
                                            PWxList=PWxList(1:length(PWxList)-1);
                                            PWyList=PWyList(1:length(PWyList)-1);
                                            signList=signList(1:length(signList)-1);
                                        end

                                    else

                                       % Elimnate the pinwheels at the borders of the ROI    
                                       count = count-1;
                                       PWxList=PWxList(1:length(PWxList)-1);
                                       PWyList=PWyList(1:length(PWyList)-1);
                                       signList=signList(1:length(signList)-1);

                                    end

                                     % Handover the variable to avoid
                                     % double counting of pinwheels from
                                     % previous loop to the next iteration
                                     doppeltestx=qx;
                                     doppeltesty=qy;
                                     doppeldet=sign(det);
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
     
    contours.cReal = cReal;
    contours.cRealIndex = cRealIndex;
    contours.cRealElements = cRealElements;

    contours.cImag = cImag;
    contours.cImagIndex = cImagIndex;
    contours.cImagElements = cImagElements;

    % Plot the pinwheel positions
    if verboseMode
        figure(fig2);
        hold on;
        for ii  = 1:length(PWxList)
            if signList(ii) > 0
                plot(PWxList, PWyList,'go', 'MarkerFaceColor', 'w');
                text(PWxList(ii), PWyList(ii),num2str(ii), 'fontsize', 12, 'color',[0 0.7 0]);
            else
                plot(PWxList, PWyList,'ko', 'MarkerFaceColor', 'w');
                text(PWxList(ii), PWyList(ii),num2str(ii), 'fontsize', 12, 'color',[0 0 0.7]);
                
            end
        end
        hold off;
        
        
    end
end


%-------------------------------------------------------------------------%
%------------------ USED FUNCTIONS ---------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [answer,m,n] = dotheycross(ax,ay,bx,by,cx,cy,dx,dy)

    m=((cx-ax)*(by-ay)-(cy-ay)*(bx-ax))/((cx-dx)*(by-ay)-(cy-dy)*(bx-ax));
    n=((cx-ax)*(cy-dy)-(cy-ay)*(cx-dx))/((bx-ax)*(cy-dy)-(by-ay)*(cx-dx));

    if ((0 <= m) && (m <= 1) && (0<=n) && (n<=1)) 
        answer = true;
    else
        answer = false;
    end

end

%-------------------------------------------------------------------------%
function [c cIndex cElements cLength] = estimateContour(map)

    c=contourc(map,[0,0]);
    cIndex = find(c(1,:) == 0); %%% contourlines that 
    cElements = c(2,cIndex); %%% Number of  (x,y) vertices in the contour line
    cLength = length(cIndex);
end


function [an,th0,ph0,sig]=pw_parameter(f)

    %Bsp [a b c d]=pw_parameter([1,2,3,1]')
    %gibt einen Vektor:
    %para(1)=an  Anisotropie d.Pws
    %para(2)=th  ueberepraesentierte Orientierung (-Pi/2,Pi/2)
    %para(3)=phi in(0,2 Pi)
    %para(4)=chirality +/- 1
    %
    %Berechnet Drehwinkel th und ph0 aus denen das PW {{1,2},{3,4}}
    %aus dem Referenz PW {{1-anisotropy,0},{0,+-1}}
    %hervorgeht. Siehe 03/30Apr_3.nb
    a=f(1,:); b=f(2,:); c=f(3,:); d=f(4,:);
    
    al  = c^2+d^2;
    ga  = a^2+b^2;
    be  = -(a*c+b*d);
    
    an  = 1-sqrt( (al-sqrt(4*be^2+(al-ga)^2)+ga) / (al+sqrt(4*be^2+(al-ga)^2)+ga));
    th0 = -atan( 2*be/(ga-al+sqrt(4*be^2+(al-ga)^2)) );
    
    M   = [a,b;c,d];
    
    vec = inv(M)*[2*be^2;be*(ga-al-sqrt(4*be^2+(al-ga)^2))];
    ph0 = angle(vec(1)+1i*vec(2));
    
    sig = sign(det(M));
   
    % rho=atan(c/a);
   % sigma=-atan(b/d)-rho;
   % amp=a/cos(rho);
   % alp=-c/(amp*sin(rho+sigma));
end