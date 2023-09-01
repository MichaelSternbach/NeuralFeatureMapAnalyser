function [z,FreeParamter] = makeOPM(MapType,X,Y,num_hc,FreeParameter)
    if isLRI(MapType)
        if size(char(MapType),2)>3
            num_amp = getNumAmp(MapType);
        else
            num_amp = 18;
        end
        if nargin == 4
            [z, FreeParamter] = makeLRI_OPM(X,Y,num_hc,num_amp);
        else
            [z, FreeParamter] = makeLRI_OPM(X,Y,num_hc,num_amp,FreeParameter);
        end
    elseif isRWM_Sym(MapType)
        if size(char(MapType),2)>7
            p_RWM = get_p_RWM_Sym(MapType);
        else
            p_RWM = 0.7;
        end
        
        if nargin == 4
            [z, FreeParamter] = makeRWM_OPM_Sym(X,Y,num_hc,p_RWM);
        else
            [z, FreeParamter] = makeRWM_OPM_Sym(X,Y,num_hc,p_RWM,FreeParameter);
        end
        
    elseif isRWM(MapType)
        if size(char(MapType),2)>3
            p_RWM = get_p_RWM(MapType);
        else
            p_RWM = 0.7;
        end
        
        if nargin == 4
            [z, FreeParamter] = makeRWM_OPM(X,Y,num_hc,p_RWM);
        else
            [z, FreeParamter] = makeRWM_OPM(X,Y,num_hc,p_RWM,FreeParameter);
        end
    else
        error('Error: MapType not known!')
    end
end

function p_RWM = get_p_RWM_Sym(MapType)
    character = char(MapType);
    p_RWM = str2double(string(character(1,8:size(character,2))));
end

function RWM = isRWM_Sym(MapType) 
    character = char(lower(MapType));
    RWM = (string(character(1:7)) == "rwm_sym");
end

function p_RWM = get_p_RWM(MapType)
    character = char(MapType);
    p_RWM = str2double(string(character(1,4:size(character,2))));
end
    
function [z, FreeParamter] = makeRWM_OPM_Sym(X,Y,num_hc,p_RWM,FreeParamter)
    

    %% randomly set free Parameters if not given
    if nargin ==4
        theta_ON = rand*pi;
        FreeParamter.theta_ON = theta_ON;
        Conj = rand>0.5;
        FreeParamter.Conj = Conj;
    elseif nargin == 5
        theta_ON = FreeParamter.theta_ON;
        Conj = FreeParamter.Conj;
    end

    % when the angle difference between ON and OFF is 5 deg,
    % 10 * r_ON equals roughly one hypercolumn
    delta_d = 0;
    delta_theta = 5*pi/180;

    % distances and angle parameters for ON and OFF grids
    
    d_ON = 128;

    theta_OFF = theta_ON + delta_theta;
    d_OFF = d_ON*(1+delta_d);

    % phases of Fourier modes
    ph0 = ( exp(1i*(theta_ON+theta_OFF)) .* ( d_ON.*exp(1i*theta_OFF) + d_OFF.* exp(1i*theta_ON) ) ) ./ ( d_ON.*exp(1i*theta_ON) + d_OFF.*exp(1i*theta_OFF));

    %kc = 4*pi/(sqrt(3)*(d_OFF*d_ON))*sqrt(sum(d_OFF^2)+sum(d_ON^2)-d_OFF*d_ON*cos(theta_OFF-theta_ON))


    ph = zeros(1,6);
    ph(1) = -1/2 * (1 + 1i*sqrt(3)) * ph0;
    ph(2) = -1/2 * (1 + 1i*sqrt(3)) * ph0;
    ph(3) = -1/2 * (1 - 1i*sqrt(3)) * ph0;
    ph(4) = -1/2 * (1 - 1i*sqrt(3)) * ph0;
    ph(5) =  ph0;
    ph(6) =  ph0;

    % Location of the modes (angle in critical circle around origin)
    er_ON = pi/d_ON          * ([cos(theta_ON),sin(theta_ON)]);
    ep_ON = pi/(sqrt(3)*d_ON)* ([-sin(theta_ON),cos(theta_ON)]);

    er_OFF = pi/d_OFF          * ([cos(theta_OFF),sin(theta_OFF)]);
    ep_OFF = pi/(sqrt(3)*d_OFF)* ([-sin(theta_OFF),cos(theta_OFF)]);

    k = zeros(2,6);
    k(:,1) =  2 *((er_ON - er_OFF) + (ep_ON - ep_OFF));
    k(:,2) = -2 *((er_ON - er_OFF) + (ep_ON - ep_OFF));
    k(:,3) =  2 *((er_ON - er_OFF) - (ep_ON - ep_OFF));
    k(:,4) = -2 *((er_ON - er_OFF) - (ep_ON - ep_OFF));
    k(:,5) =  4 *(ep_ON - ep_OFF);
    k(:,6) = -4 *(ep_ON - ep_OFF);

    %% scale map
    map_dim = size(X,1);
    kc = sqrt(sum(k(:,1).^2));
    %pixels_per_hc =map_dim/num_hc;
    scale = 1/kc*2*pi/map_dim*num_hc;
    k=k.*scale;

    %i = reshape(1:6,[1,1,6]);
    ph = reshape(ph,[1,1,6]);

    %k=k./100;
    kx = reshape(k(1,:),[1,1,6]);
    ky = reshape(k(2,:),[1,1,6]);

    z = sum(exp(1j*(X.*kx+Y.*ky)).*ph,3);

    %% normalise
    z = z - mean(z(:));
    z = z/sqrt(mean(abs(z(:)).^2));
    
    %% add Gaussian random field
    z = z*p_RWM + makeGaussianRandomField(map_dim,num_hc).*(1-p_RWM);
    
    if Conj
        z = conj(z);
    end


end

function [z, FreeParamter] = makeRWM_OPM(X,Y,num_hc,p_RWM,FreeParamter)
    

    %% randomly set free Parameters if not given
    if nargin ==4
        theta_ON = rand*pi;
        FreeParamter.theta_ON = theta_ON;
    elseif nargin == 5
        theta_ON = FreeParamter.theta_ON;
    end

    % when the angle difference between ON and OFF is 5 deg,
    % 10 * r_ON equals roughly one hypercolumn
    delta_d = 0;
    delta_theta = 5*pi/180;

    % distances and angle parameters for ON and OFF grids
    
    d_ON = 128;

    theta_OFF = theta_ON + delta_theta;
    d_OFF = d_ON*(1+delta_d);

    % phases of Fourier modes
    ph0 = ( exp(1i*(theta_ON+theta_OFF)) .* ( d_ON.*exp(1i*theta_OFF) + d_OFF.* exp(1i*theta_ON) ) ) ./ ( d_ON.*exp(1i*theta_ON) + d_OFF.*exp(1i*theta_OFF));

    %kc = 4*pi/(sqrt(3)*(d_OFF*d_ON))*sqrt(sum(d_OFF^2)+sum(d_ON^2)-d_OFF*d_ON*cos(theta_OFF-theta_ON))


    ph = zeros(1,6);
    ph(1) = -1/2 * (1 + 1i*sqrt(3)) * ph0;
    ph(2) = -1/2 * (1 + 1i*sqrt(3)) * ph0;
    ph(3) = -1/2 * (1 - 1i*sqrt(3)) * ph0;
    ph(4) = -1/2 * (1 - 1i*sqrt(3)) * ph0;
    ph(5) =  ph0;
    ph(6) =  ph0;

    % Location of the modes (angle in critical circle around origin)
    er_ON = pi/d_ON          * ([cos(theta_ON),sin(theta_ON)]);
    ep_ON = pi/(sqrt(3)*d_ON)* ([-sin(theta_ON),cos(theta_ON)]);

    er_OFF = pi/d_OFF          * ([cos(theta_OFF),sin(theta_OFF)]);
    ep_OFF = pi/(sqrt(3)*d_OFF)* ([-sin(theta_OFF),cos(theta_OFF)]);

    k = zeros(2,6);
    k(:,1) =  2 *((er_ON - er_OFF) + (ep_ON - ep_OFF));
    k(:,2) = -2 *((er_ON - er_OFF) + (ep_ON - ep_OFF));
    k(:,3) =  2 *((er_ON - er_OFF) - (ep_ON - ep_OFF));
    k(:,4) = -2 *((er_ON - er_OFF) - (ep_ON - ep_OFF));
    k(:,5) =  4 *(ep_ON - ep_OFF);
    k(:,6) = -4 *(ep_ON - ep_OFF);

    %% scale map
    map_dim = size(X,1);
    kc = sqrt(sum(k(:,1).^2));
    %pixels_per_hc =map_dim/num_hc;
    scale = 1/kc*2*pi/map_dim*num_hc;
    k=k.*scale;

    %i = reshape(1:6,[1,1,6]);
    ph = reshape(ph,[1,1,6]);

    %k=k./100;
    kx = reshape(k(1,:),[1,1,6]);
    ky = reshape(k(2,:),[1,1,6]);

    z = sum(exp(1j*(X.*kx+Y.*ky)).*ph,3);

    %% normalise
    z = z - mean(z(:));
    z = z/sqrt(mean(abs(z(:)).^2));
    
    %% add Gaussian random field
    z = z*p_RWM + makeGaussianRandomField(map_dim,num_hc).*(1-p_RWM);
    


end


function RWM = isRWM(MapType) 
    character = char(lower(MapType));
    RWM = (string(character(1:3)) == "rwm");
end

function num_amp = getNumAmp(MapType)
    character = char(MapType);
    num_amp = str2num(string(character(1,4:size(character,2))));
end

function LRI = isLRI(MapType) 
    character = char(lower(MapType));
    LRI = (string(character(1:3)) == "lri");
end

function [z, FreeParameter] = makeLRI_OPM(X,Y,num_hc,num_amp,FreeParameter)
    
    %% Scale map
    map_dim = size(X,1);
    pixels_per_hc =map_dim/num_hc;
    kc=2*pi/pixels_per_hc;

    %% randomly set free Parameters if not given
    if nargin ==4
        li = (-1).^ceil(rand(1,1,num_amp)-0.5);
        phii= 2*pi*rand(1,1,num_amp);
        
        FreeParameter.li = li;
        FreeParameter.phii = phii;
    elseif nargin == 5
        li = FreeParameter.li;
        phii = FreeParameter.phii;
    end
    
    
    %% superposition planewaves
    i = reshape(1:num_amp,[1,1,num_amp]);
    z = sum(exp(1j*((kc*li.*(X.*cos(o_n(i,num_amp))+Y.*sin(o_n(i,num_amp))))+phii)),3);
    
    %% normalise
    z = z - mean(z(:));
    z = z/sqrt(mean(abs(z(:)).^2));
    
end

function o = o_n(i,num_amp)
    o = i.*pi/num_amp;
end