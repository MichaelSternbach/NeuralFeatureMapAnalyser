function [local_pw_dens, local_pw_dens_plus, local_pw_dens_minus] = put_gaussians(sizex,sizey, PWxList, PWyList, Lambda, sig, roi_matrix, sign_list)
%
% FUNCTION IS CALLED BY ANALYZE_SINGLE_MAP.M
%
% INPUT ARGUMENTS
%
% sizex      ... x-size of data matrix  
% sizey      ... y-size of data matrix   
% Lambda     ... average wavelength in pixels
% measure    ... how many pixels correspond to one millimeter
% sig        ... width of the smoothing kernel in units of Lambda
% roi_matrix ... region of interest, matrix being 1 inside the ROI and 0
%                otherwise
% sign_list  ... list of signs of the pinwheels
% choose 1 for fully automated method according to Kaschube et al.,Science, 2010, Supp. Inf.
% choose 0.25 for semi-automated method according to Kaschube et al.,Science, 2010, Supp. Inf.
% choose 0.01 for pw-position estimation according to Kaschube et al.,Science, 2010, Supp. Inf.
%
%
% OUTPUT ARGUMENTS
%
% local_pw_dens      ... local pinwheel density with Gaussians put on each pinwheel
%                        position
% local_pw_dens_plus ... local pinwheel density for positively charged pinwheels with Gaussians put on each pinwheel
%                        position
% local_pw_dens_minus... local pinwheel density for negatively charged pinwheels with Gaussians put on each pinwheel
%                        position
%
%
    local_pw_dens = zeros(sizex,sizey);
    local_pw_dens_plus = zeros(sizex,sizey);
    local_pw_dens_minus = zeros(sizex,sizey);
    [X Y] = meshgrid(1:sizey,1:sizex);
    
    sig = sig*Lambda;

    for ii = 1:length(PWxList)
        local_pw_dens = local_pw_dens + 1/2/pi/sig^2*exp(-1/2/sig^2*((X-PWxList(ii)).^2+(Y-PWyList(ii)).^2));    
        if nargout > 1
            if nargin < 8
                error('Please pass on PW-signList to function if more than one output arg is specified');
            else  
                if sign_list(ii) > 0
                    local_pw_dens_plus = local_pw_dens_plus + 1/2/pi/sig^2*exp(-1/2/sig^2*((X-PWxList(ii)).^2+(Y-PWyList(ii)).^2));    
                else
                    local_pw_dens_minus = local_pw_dens_minus + 1/2/pi/sig^2*exp(-1/2/sig^2*((X-PWxList(ii)).^2+(Y-PWyList(ii)).^2));    
                end
            end
        end
    end
    
    %%%% Rescale with the overlap    
    ext =round(2.^(ceil(log(max([sizex,sizey]))/log(2))));   %(This will be the size of the filters/arrays)
    a=[sizex,sizey];
    
    roi_padd = zeros(ext);
    local_pw_dens_padd = zeros(ext);

    roi_padd(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)))=roi_matrix;
    local_pw_dens_padd(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)))=local_pw_dens;
    
    [xx,yy] = meshgrid(1:ext,1:ext);   % create coordinate system
    xx=xx-ext/2;
    yy=yy-ext/2;
    
    filt= 1/2/pi/sig^2 * exp(-((xx.^2+yy.^2)/2/sig^2));
    overlap = fftshift(ifft2(fft2(fftshift(roi_padd)).*fft2(fftshift(filt))));                
    local_pw_dens_padd(roi_padd == 1) = local_pw_dens_padd(roi_padd == 1)./overlap(roi_padd == 1);
    
    local_pw_dens=local_pw_dens_padd(round((ext-a(1))/2+1):round((ext-a(1))/2+a(1)),round((ext-a(2))/2+1):round((ext-a(2))/2+a(2)));

    
end %%%% END OF FUNCTION PUT_GAUSSIANS.M
