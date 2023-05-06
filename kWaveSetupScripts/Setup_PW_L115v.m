%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --kWave Script for Plane Wave Transmission (L11-5v transducer)
% --Author: Gayathri M and Mahesh Raveendranatha Panicker
% --Date: 04-05-2023
% --Center for Computational Imaging, Indian Institute of Technology Palakkad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;

filepath = pwd;
addpath(genpath(pwd));

dataset_name = 'Specular_Reflector_20';

GPU             = false;   % default running on CPU
if (gpuDeviceCount)
    GPU         = true;   % set to true to simulate using GPU, to false to simulate in CPU
end

USE_CPU_C_CODE  = false;  % set to true to simulate using CPU cores in the absence of CUDA
SimPlot         = false; % set to true to view wave propagation. Will be active only if USE_CPU_CORES == false
RecordMovie     = false; % set to true to record the movie of simulation, else set to false
DO_BEAMFORMING  = true;  % set to true to perform beamforming with the generated dataset, else set to false

%% Define L11-5V transducer
Trans.name = 'L11-5v';                % Transducer name
Trans.tx_N_elements = 128;            % Number of transducer elements
Trans.tx_element_width = 270e-6;      % Width of a single transducer element [m]
Trans.tx_kerf = 30e-6;                % Kerf of the transducer [m]
Trans.tx_pitch = 300e-6;              % Pitch of the transducer [m]
Trans.tx_elevation_height = 5e-3;     % Elevation height of the transducer [m]
Trans.tx_elevation_focus = 18e-3;     % Elevation focus of the transducer [m]
Trans.tx_Fc = 7.6e6;                  % Transmit frequencies [Hz]
Trans.tx_Fmax = 1.2*Trans.tx_Fc;      % Maximum possible transducer center frequency [Hz]
Trans.tx_c0 = 1540;                   % SoS of the transducer [m/s]

Trans.tx_width = Trans.tx_N_elements * Trans.tx_pitch;% Total width of the transducer [m]
Trans.na = 7;                          % Number of plane wave steering angles
Trans.angles = linspace(-18,18,Trans.na);     % Steering angle of the transmit signal between [-18, +18] [degrees]
probe_geometry = linspace(-Trans.tx_width/2, Trans.tx_width/2, Trans.tx_N_elements)';

%% DEFINE THE K-WAVE GRID

endDepth = 0.7*3e-2;     % Maximum depth in diagonal dir [m](1/sqrt(2) = 0.7)
points_per_wavelength = 4;

% Set the size of the perfectly matched layer (PML)
PML_Z_size = 20;
PML_X_size = 20;

c_min = Trans.tx_c0;
dx = c_min/(points_per_wavelength*Trans.tx_Fmax);% grid spacing along lateral direction [m]
dz = dx;   % grid spacing along axial direction [m]
               
Nz = ceil(endDepth/dz); % grid points along z
Nx = ceil(Trans.tx_width/dx); % grid points along x

% Rounding grid points to nearest power to reduce computations
Nz = 2^nextpow2(Nz)- 2*PML_Z_size;
Nx = 2^nextpow2(Nx)- 2*PML_X_size;

% Create the k-space grid
kgrid = kWaveGrid(Nz,dz, Nx, dx);

% Create the time array
t_end = (Nz * dz) * 2.2 / Trans.tx_c0;   % [s]
CFL = 0.05;
kgrid.t_array = makeTime(kgrid, Trans.tx_c0, CFL, t_end);

%% DEFINE THE MEDIUM PARAMETERS

% Mimicking a specular reflector inclined at an angle 20 degree
lineStartZ = round(Nz/2);    %Start X grid point of the line
lineStartX = round(Nx/2);    %Start Y point of the line
lineLen = 300;               %Length of the line
lineThickness = 8;           %Thickness of the line
tilt = -20;                  %Tilt of the line

% Phantom function
medParam  = specularReflector (Nz, Nx, lineStartZ, lineStartX, lineLen, lineThickness, tilt);

medium.sound_speed            = medParam.sound_speed;
medium.density                = medParam.density;
medium.alpha_coeff            = medParam.alpha_coeff;
medium.alpha_power            = 1.5;
medium.BonA                   = 6;

ax_vec  = kgrid.x_vec - min(kgrid.x_vec);
lat_vec = kgrid.y_vec - min(kgrid.y_vec);

figure,imagesc(ax_vec.*1000, lat_vec.*1000,medium.sound_speed); % Plotting the speed of sound map
title('SoS map');
colorbar;
figure,imagesc(ax_vec.*1000, lat_vec.*1000,medium.density); % Plotting the density map
title('Density map');
colorbar;

%% DEFINE THE INPUT SIGNAL
source_strength     = 4e6;      % [Pa]
source_cycles       = 3;        % number of tone burst cycles
source_freq         = Trans.tx_Fc;   
source.p_mask   = zeros(Nz, Nx);
source_width    = round(Trans.tx_element_width/dz);
source_kerf     = round(Trans.tx_kerf/dz);
element_height  = round(Trans.tx_elevation_height/dz);
rotation        = 0;

% Create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

% Creating each transducer element by combining grid points
for N_c = 1: Trans.tx_N_elements
    karray.addRectElement([kgrid.x_vec(1), probe_geometry(N_c, 1)],Trans.tx_kerf, Trans.tx_element_width, rotation);
end

% Defining source and sensor arrays
karray.setArrayPosition([0, 0], 0);
source.p_mask = karray.getArrayBinaryMask(kgrid);

display_mask=medium.sound_speed/max(medium.sound_speed(:));
display_mask=histeq(display_mask);
display_mask(display_mask<0.7)=0;
display_mask(display_mask>0.7)=1;

sensor.mask = zeros(Nz, Nx);
sensor.mask = source.p_mask;

%% RUN THE SIMULATION

acq = 'PW';
for txIdx = 1: length(Trans.angles)
    
    disp(txIdx);
    
    % Calculate the steering delay
    if(Trans.angles(txIdx)>0)
        Trans.delay(:,txIdx) = (probe_geometry(Trans.tx_N_elements,1)-probe_geometry(1:Trans.tx_N_elements,1))*sind(Trans.angles(txIdx))/(Trans.tx_c0*kgrid.dt);
    else
        Trans.delay(:,txIdx) = (probe_geometry(1,1)-probe_geometry(1:Trans.tx_N_elements,1))*sind(((Trans.angles(txIdx))))/(Trans.tx_c0*kgrid.dt);
    end
    
    % Note that the source strength is not scaled by acoustic impedance as a pressure source is used (Z=P/V)
    source_sig = source_strength*toneBurst(1/kgrid.dt, source_freq, source_cycles,'SignalOffset',Trans.delay(:,txIdx));
    source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
    
    if(GPU)
        DATA_CAST       = 'gpuArray-single';
        input_args = {'PMLSize', [PML_Z_size, PML_X_size], 'PlotPML', false, ...
            'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
            'DisplayMask', display_mask, 'DataCast', DATA_CAST, 'DataRecast', true,...
            'RecordMovie', RecordMovie, 'MovieProfile', 'MPEG-4'};
        if(USE_CPU_C_CODE == true)
            sensor_data_grid_points = kspaceFirstOrder2DC(kgrid, medium, source, sensor, input_args{:});
        elseif(SimPlot == true)
            sensor_data_grid_points = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        else
            sensor_data_grid_points = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
        end
    else
        DATA_CAST       = 'single';
        if(SimPlot == false)
            input_args = {'PMLSize', [PML_Z_size, PML_X_size], 'PMLInside', false, 'DataCast', DATA_CAST,  'DataRecast', true,'PlotSim', false};
        else
            input_args = {'PMLSize', [PML_Z_size, PML_X_size], 'PlotPML', false, ...
            'PMLInside', false, 'PlotScale', [-1, 1]*source_strength, ...
            'DisplayMask', display_mask, 'DataCast', DATA_CAST, 'DataRecast', true,...
            'RecordMovie', RecordMovie, 'MovieProfile', 'MPEG-4'};
        end
        sensor_data_grid_points = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    end
    sensor_data(:,:,txIdx) = karray.combineSensorData(kgrid, sensor_data_grid_points); % Combine the grid points to get in terms of physical transducer elements
    
end

RcvData = permute(sensor_data, [2 1 3]);

% Saving the dataset
dataset = dataExtractFunc(dataset_name,acq,Trans,probe_geometry,RcvData,kgrid,medium);

%%
if(DO_BEAMFORMING)
    
    % DAS PW Beamforming
    Fnum = 1.5; % Receive F# for beamforming

    acq = 'PW';

    rawData=double(dataset.rawData);
    Fs = dataset.Fs*1e6;
    timeVector = (0:(size(rawData,1)-1)).'/(Fs);
    
    z_axis = 0.5*1540*timeVector;
    [x_grid,z_grid] = meshgrid(dataset.probe_geometry,z_axis);
    
    x_lim = [min(x_grid(:)) max(x_grid(:))]*1e3;
    z_lim = [min(z_grid(:)) max(z_grid(:))]*1e3;
    
    ang = dataset.angles;

    for txIdx = 1: length(ang)
        beamformedDataPW (:,:,txIdx) = DAS_RF(acq,squeeze(rawData(:, :, txIdx)), timeVector, x_grid, z_grid, dataset.probe_geometry, -ang(txIdx), dataset.Trans.tx_c0*ones(size(x_grid)), Fnum);
    end
    
    beamformedDataDAS = sum(beamformedDataPW(:,:,:), 3);
    envelopeDAS = abs(hilbert(beamformedDataDAS));
    beamformedDataDASImage = (envelopeDAS(:,:)./max(max(envelopeDAS(:,:))));
    figure,imagesc(dataset.probe_geometry.*1000,z_axis.*1000,20*log10(beamformedDataDASImage));
    colormap(gray);
    colorbar;
    vrange = [-60 0];
    caxis(vrange);
    xlabel('x [mm]');
    ylabel('z [mm]');
    title('DAS-PW');
    set(gca,'fontsize',20);
    axis([x_lim z_lim]);
    
    fsavepath = strcat(filepath,filesep,'kWaveImages',filesep,dataset_name);
    if ~exist(fsavepath, 'dir')
        mkdir(fsavepath);
    end
    
    saveas(gcf,[fsavepath,filesep,dataset_name,'_',acq,'_',dataset.timeStamp,'.fig']);
    
end

