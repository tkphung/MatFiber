function [XYAI,varargout] = sfvectors(SFIMG,varargin)
%sfvectors(SFIMG,varargin): Segmentation of stress fibers developed
%from MatFiber (Fomovsky & Holmes 2010). This generates a fiber direction
%in per set square pixel window based on a gradient-intensity method.
% 
% [XYAI,varargout] = stressfiberseg(SFIMG,varargin)
%   INPUTS:
%       SFIMG- StresS Fiber Image file (phalloidin stained, 20x mag)
%       varargin: 'visualize' followed by true or false (default)
%                 'pixres' followed by microns/pixel resolution
%                          (0.29 um/pixel default- 20x mag)
%                 'threshold' followed by pixel intensity to threshold
%                 'window' followed by micron window size (3 um default)
%                 'borderoverlap' followed by border logical mask
%                 'borderthresh' followed by % acceptable border overlap
%   OUTPUTS:
%       XYAI- x position, y position, fiber angle, mean window pixel
%             intensity
%       varargout: bcount- % pixels in non- basal cell body regions
% 
% References used:
% This code is based off of MatFiber (Fomovsky & Holmes 2010)
%
% Created by Thien-Khoi N. Phung (April 17, 2020)

% Deal with VARARGIN
visual_flag  = false;
pr           = 0.29;  % default pixel resolution (um/pixel)
window       = 3;     % window for vectors (um)
thresh_flag  = false;
border_flag  = false;
bthresh_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'visualize'
                visual_flag = varargin{jz+1};
            case 'pixres'
                pr = varargin{jz+1};
            case 'threshold'
                thresh_flag = true;
                pthresh     = varargin{jz+1};
            case 'window'
                window = varargin{jz+1};
            case 'borderoverlap'
                border_flag = true;
                border      = varargin{jz+1};
            case 'borderthresh'
                bthresh_flag = true;
                bthresh      = varargin{jz+1};
            otherwise
                error('SF ERROR: Check your varargins.')
        end
    end
end % varargin

% Convert image to double
    dImage = double(SFIMG);
    
% SUBREGION SIZE (S) determines how many pixels (S x S) you wish to group 
% into only 1 output angle.
% The code will divide the original image into N1 x N2 subregions of SxS
% pixels each, and calculate 1 orientation value for each subregion. This
% orientation will correspond to the angle perpendicular to the strongest 
% mean intensity gradient direction within that subregion.
% Determine subregion size based on pixres (to normalize across low and
% high magnification images)
    S = round(window/pr); % pixels for a # um window (rounded)
    
% Initial computations related to the image and subregions
% Because gradient mask (next step) is 9x9, convolution result will be
% erroneous in pixels near image border. Subregions are shifted 5 pixels
% away from the image borders to avoid this error.
    D  = size(dImage);
    N1 = floor((D(1)-10)/S);           % Number of subregions vertically;
    N2 = floor((D(2)-10)/S);           % Number of subregions horizontally;

% Gradient masks computation
    for ii = -4:4
        for jj = -4:4
            MX(ii+5,jj+5) = 2*ii/9*exp(-(ii^2+jj^2)/4);
            MY(ii+5,jj+5) = 2*jj/9*exp(-(ii^2+jj^2)/4);
        end
    end
    
% Convolve masks with the image
    GX = conv2(dImage,MX,'same');
    GY = conv2(dImage,MY,'same');

% Edge image and gradient direction
    E = GX.^2+GY.^2;
    phi = 180/pi*atan2(GY, GX);

% Determine local orientation in every subregion
    bins       = 0:1:179;
    FiberAngle = nan([N1*N2 1]);
    FiberPosX  = FiberAngle;
    FiberPosY  = FiberAngle;
    mpi        = FiberAngle;
    bcount     = FiberAngle;
    count = 0;
    for qq = 1:1:N1
        for rr = 1:1:N2
            % Subregion window to analyze
            lx = 5 + S*(rr - 1) + 1;
            ux = 5 + S*rr;
            ly = 5 + S*(qq - 1) + 1;
            uy = 5 + S*qq;

            % AthetaW is a sum of alignment values between pixel gradient 
            % directions and each direction from 'bins', weighted by pixel
            % gradient magnitudes.
            AthetaW = zeros(size(bins));
            for nn = 1:1:numel(bins)
                C = E(ly:uy,lx:ux).*exp(2*cos(2*pi/180*(bins(nn) - phi(ly:uy,lx:ux))))/exp(2);
                AthetaW(nn) = sum(sum(C));
            end
            
             % Sets 'Ang' as the entry from 'bins' corresponding to the maximum 'AthetaW' value.
            [~, Ang] = max(AthetaW);

            % Adjust between -90 and 90 degrees
            if Ang > 90
                Ang = Ang - 180;
            end

            % Save the fiber direction and centroid of subregion
            count = count + 1;
            FiberAngle(count) = Ang;
            FiberPosY(count) =  round(5 + S*qq - S/2);
            FiberPosX(count) =  round(5 + S*rr - S/2);
            
            % mean pixel intensity in subregion
            pxls       = SFIMG(ly:uy,lx:ux);
            mpi(count) = mean(pxls(:));
            
            if border_flag
                % Count pixels in border region (% of total window pixels)
                bpxls         = border(ly:uy,lx:ux);
                bcount(count) = sum(bpxls(:))./numel(bpxls);
            end
            
            clear AthetaW Ang;
        end
    end

% ALl of the SF
SFXYAIog = [FiberPosX(:) FiberPosY(:) FiberAngle(:) mpi(:)];
OGbcount = bcount;

% (OPTION) Threshold based on basal cell boundaries
if thresh_flag    
    % Which vectors to include?
    keepvecs = mpi>=pthresh;
    
    FiberPosX(~keepvecs)  = [];
    FiberPosY(~keepvecs)  = [];
    FiberAngle(~keepvecs) = [];
    mpi(~keepvecs)        = [];
   
    if border_flag
        bcount(~keepvecs)     = [];
    end
end

% (OPTION) Threshold based on basal cell boundaries overlap percentage
if bthresh_flag    
    % Which vectors to include?
    keepvecs = bcount<=bthresh;
    
    FiberPosX(~keepvecs)  = [];
    FiberPosY(~keepvecs)  = [];
    FiberAngle(~keepvecs) = [];
    mpi(~keepvecs)        = [];
    bcount(~keepvecs)     = [];
end

% OUTPUT: Position-Angle-Mean Intensity
XYAI = [FiberPosX(:) FiberPosY(:) FiberAngle(:) mpi(:)];

% (OPTION) Percentage of pixels in "border"
if border_flag
    varargout{1} = OGbcount;
    varargout{2} = SFXYAIog;
end

% (OPTION) Visualize fiber vectors
if visual_flag
    % Calculate the Fiber Vectors
    FiberOrientationX = cos(FiberAngle*pi/180);
    FiberOrientationY = sin(FiberAngle*pi/180);
    
    figure('WindowStyle','docked','NumberTitle','off','name','SF Vectors')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(imadjust(SFIMG),'Parent',axt);
    hold on;
    quiver(FiberPosX,FiberPosY,FiberOrientationX,-FiberOrientationY,0.3,'y','LineWidth',1.5,...
        'ShowArrowHead','off');
    axis equal tight
end % visualize

end % sfvectors()