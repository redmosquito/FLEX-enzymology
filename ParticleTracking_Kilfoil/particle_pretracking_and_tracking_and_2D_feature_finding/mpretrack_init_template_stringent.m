%% mpretrack_init_template 
% Feature Finding mpretrack_init Template
% Drew Drabek, Blacklow Lab, Harvard Medical School, BCMP
% Based on code from http://people.umass.edu/kilfoil/tools.php
% Maria Kilfoil at UMass Amherst

%% Set Parameters for pretracking
clear;clc;  %initialize

% navigate to the folder with your images
basepath = [cd,'/'];  %set basepath

% Run the script mpretrack_renameFOV.m

% Define your parameters
featsize = 7;   %radius of reatures in pixels, The size of the feature you want to find.
barint = 5;   %The minimum intensity you want to accept.
barrg = 40;      %The maximum Rg squared you want to accept.
barcc = .25;    %The maximum eccentricity you want to accept.
IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
fovn = 0;       %ID# for the series of images (typically one field of view)
frame = 0;      %ID# for the frame in the series
Imin=1000;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
masscut=0;      %(optional) the masscut parameter for feature2D to remove false positives before rifining the position to speed up the code.
field=2;        %(optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
                %    image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
                %    set to 2 for progressive scan cameras. Defaults to 2.

%% Run pretracking function and save output
%Run the algorithm with the defined parameters
tic
[M2, MT] = mpretrack_init(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame, Imin, masscut, field);
toc

%save the output
mkdir('mpretrack_init')
savefig(fullfile(basepath,'mpretrack_init','fov0_pretrackfig'));
close(gcf);
save(fullfile(basepath,'mpretrack_init','fov0_pretrackparams'));

%% Feature finding loop function for all fields of view
%calculate the number of frames
numframes=length(dir(fullfile(['fov',num2str(fovn)],'*.tif')))-1;
%numframes=5;
%pretrack all the frames specified in numframes
tic
mpretrack(basepath, fovn, featsize, barint, barrg, barcc, IdivRg, numframes, Imin, masscut, field);
toc

%%