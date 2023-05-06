%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines a specular reflector in a homogeneous medium
% Inputs:
%         Nx - grid points in x direction
%         Ny - grid points in y direction
%         lineStartX - Starting grid point of the reflector (lateral)
%         lineStartZ - Starting grid point of the reflector (axial)
%         lineLen - Length of the reflector (grid points)
%         lineThickness - Thickness of the reflector (grid points)
%         tilt - Orientation of the reflector (degree)
% Author: Gayathri Malamal
% Date: 04-05-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function medParam = specularReflector(Nz, Nx,lineStartZ, lineStartX, lineLen, lineThickness, tilt)
%==========================================================================
% Define the medium properties
% =========================================================================
% Tissue 
cTissue                 = 1540;   % Avg speed of sound in tissue [m/s]
rhoTissue               = 1000;   % Density [kg/m^3]
alphaTissue             = 0.0022; % Absorption [dB/(MHz cm)]

% Reflector
cReflector              = 3000;   % ~Speed of sound in  [m/s]
rhoReflector            = 1900;   % ~Density in [kg/m^3] of bone
alphaReflector          = 3.5;    % ~Absorption [dB/(MHz cm)]

% =========================================================================
% Creating reflector 
% =========================================================================

phantom =  makeLine(Nz, Nx, [round(lineStartZ) round(lineStartX)], -tilt*pi/180, lineLen); 

for Nc=1:lineThickness
    phantom = phantom + makeLine(Nz, Nx, [round(lineStartZ)+Nc, round(lineStartX)], -tilt*pi/180, lineLen);
end

phantom = permute(phantom,[1,3,2]);

%=========================================================================
%Creating background maps of random scatterers
%=========================================================================
rng(1); % Seed initialization to ensure same scatterer distribution in every dataset
scat_dist = randn(Nz,Nx); % Random distribution of scatterers 

background_map_mean = 1;
background_map_std = 0.01;
background_map = background_map_mean + background_map_std * scat_dist;

% =========================================================================
% Defining medium speed of sound and density
% =========================================================================

medParam.sound_speed               = cTissue*ones(Nz, Nx).*background_map;
medParam.sound_speed(phantom == 1) = cReflector;
medParam.density                   = rhoTissue*ones(Nz, Nx).*background_map;
medParam.density(phantom== 1)      = rhoReflector;

medParam.alpha_coeff               = alphaTissue*ones(Nz, Nx);
medParam.alpha_coeff(phantom == 1) = alphaReflector;
end
