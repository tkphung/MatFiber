function [HETAUC,HETMET] = sfhetmetric(SFXYAIB,pr)
%sfhetmetric(SFXYAIB,pr): Takes in stress fiber vector field from
%sfvectors() and calculates alignment heterogeneity based on (Richardson &
%Holmes 2016).
% 
% [HETAUC,HETMET] = sfhetmetric(SFXYAIB,pr)
%   INPUTS:
%       SFXYAIB- Stress fiber vector field including columns for at least
%                   X coordinate
%                   Y coordinate
%                   A fiber angle
%       pr     - pixel resolution (um/pixel)
%   OUTPUTS:
%       HETAUC - heterogenity alignment under the curve (see details in
%                Richardson & Holmes 2016)
%       HETMET - histogram data used for calculation of HETAUC
%                   Column 1: distance between fibers
%                   Column 2: mean fiber-fiber alignment of sample
%                   Column 3: mean fiber-fiber alignment of randomized
%                             control
% 
% References used:
% This code is based off (Richardson & Holmes 2016)
%
% Created by Thien-Khoi N. Phung (April 17, 2020)

% Calculate heterogeneity metric
    % Calculate distance between each stress fiber and all of its nbrs
    % sfd- unique pair differences arranged in the order (2,1), (3,1), ...,
    % (m,1), (3,2), ..., (m,2), ..., (m,m–1), i.e., the lower-left triangle
    % of the m-by-m distance matrix in column order
    sfd   = pdist(SFXYAIB(:,1:2),'euclidean').*pr; % um
    
    % Keep only the pairs within a maximum distance
    maxd  = 200; % um
    keepp =  sfd<= maxd;
    sfd(~keepp) = [];
    
    % Define the indices for the kept pairs
    m            = size(SFXYAIB,1);  % Number of SF
    t            = tril(ones(m),-1); % Lower triangle
    t(t==1)      = keepp;            % Keep only the indices from keepp
    [hmid,nbrid] = find(t);          % Find hm and nbr index
    
    % Calculate the dot product for the pairs
    % Note: we 2x the angle to account for bidirectionality (-90 and 90
    % degrees are the same direction); also, we scale the dot product to be
    % between 0 and 1 (instead of -1 to 1)
    % 0 means the two vectors are perpendicular
    % 1 means the two vectors are perfectly aligned
    hmangle  = 2.*SFXYAIB(hmid,3).*(pi/180);  % radians
    nbrangle = 2.*SFXYAIB(nbrid,3).*(pi/180); % radians
    sfdot    = (1/2).*(cos(hmangle-nbrangle)+1);
    
    % CONTROL- shuffle angles and calculate dot product
    % Shuffle angles
    shuffang = SFXYAIB(randperm(size(SFXYAIB,1)),3);
    hmangle  = 2.*shuffang(hmid).*(pi/180);  % radians
    nbrangle = 2.*shuffang(nbrid).*(pi/180); % radians
    cdot     = (1/2).*(cos(hmangle-nbrangle)+1);
    
    distvals = linspace(0,maxd,50);
    HETMET   = zeros(numel(distvals)-1,3);
    for dd = 1:(numel(distvals)-1)
        
        % Which pairs are in this distance window
        pairs = sfd>=distvals(dd) & sfd<distvals(dd+1);
        
        % Average the dot products
        HETMET(dd,1) = mean(distvals(dd:(dd+1))); % distance (um)
        HETMET(dd,2) = mean(sfdot(pairs));        % mean dot product
        HETMET(dd,3) = mean(cdot(pairs));         % mean of control dot
    end
    HETMET(:,4) = repmat(mean(HETMET(:,3)),size(HETMET(:,3)));
    
    % Calculate AUC
    % Bin width
    binwidth = diff(HETMET(1:2,1));
    
    % Difference between bin heterogeneity score and the control score
    heights  = HETMET(:,2) - HETMET(:,4);
    
    % Find where the heights go negative
    negsign  = find(heights<1e-3,1);
    
    % Calculate area under the curve
    HETAUC = sum(binwidth.*heights(1:(negsign-1)));
end