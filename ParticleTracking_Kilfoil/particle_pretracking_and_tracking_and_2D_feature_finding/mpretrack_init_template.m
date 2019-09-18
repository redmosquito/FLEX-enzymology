%% mpretrack_init_template 
% Feature Finding mpretrack_init Template
% Based on code from http://people.umass.edu/kilfoil/tools.php
% Maria Kilfoil at UMass Amherst
%
% Drew Drabek, Blacklow and Loparo Labs, Harvard Medical School, BCMP

%% Set Parameters for pretracking
clear;clc;              %initialize
basepath = [cd,'\'];    % Document the working directory for the images

fovn=1;

%% Run the script mpretrack_renameFOV.m
renameFOVBF(fovn,'tif');

%% Define your parameters
featsize = 7;   %radius of reatures in pixels, The size of the feature you want to find.
barint = 2000;   %The minimum intensity you want to accept.
barrg = 100;      %The maximum Rg squared you want to accept.
barcc = 0.5;    %The maximum eccentricity you want to accept.
IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
%fovn = 0;       %ID# for the series of images (typically one field of view)
frame = 0;      %ID# for the frame in the series
Imin=50;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
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
savefig(fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_pretrackfig']));
close(gcf);
save(fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_pretrackparams']));
%% Make histograms of initial tracking
%

% plot and save histogram of particle intensities
figure; 
hist(M2(:,3),100)
    T=['fov' num2str(fovn) ' Intensity Histogram '];
    title(T);
    xlabel('Intensity');
    ylabel('count');
    savefig(fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_IntensityHistrogam.fig']));
    print('-f',fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_IntensityHistrogam.png']),'-dpng');
close(gcf);

% plot and save histogram of particle Radius of Gyration
figure; 
hist(M2(:,4),100)
    T=['fov' num2str(fovn) ' Rg Radius of Gyration Histogram '];
    title(T);
    xlabel('Radius of Gyration, Rg');
    ylabel('count');
    savefig(fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_RgHistogram.fig']));
    print('-f',fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_RgHistogram.png']),'-dpng');
close(gcf);

% plot and save histogram of particle eccentricities
figure; 
hist(M2(:,5),100)
    T=['fov' num2str(fovn) ' Eccentricity Histogram '];
    title(T);
    xlabel('Eccentricity');
    ylabel('count');
    savefig(fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_EccentricityHistogram.fig']));
    print('-f',fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_EccentricityHistogram.png']),'-dpng');
close(gcf);

% plot and save histogram of intensity over radius of gyration
figure; 
hist(M2(:,3)./M2(:,4),100)
    T=['fov' num2str(fovn) ' IdivRg Histogram '];
    title(T);
    xlabel('IdivRg');
    ylabel('count');
    savefig(fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_IdivRgHistogram.fig']));
    print('-f',fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_IdivRgHistogram.png']),'-dpng');
close(gcf);

% % save the mpretrack_init figure as a png file
% openfig(fullfile(basepath, 'mpretrack_init/', ['fov' num2str(fovn) '_pretrackfig.fig']));
%     print('-f',fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_pretrackfig.png']),'-dpng');
% close(gcf);

%% Feature finding loop function for all fields of view

% %calculate the number of frames
% numframes=length(dir(fullfile(['fov',num2str(fovn)],'*.tif')))-1;
% %numframes=100;
% %pretrack all the frames specified in numframes
% tic
% mpretrack(basepath, fovn, featsize, barint, barrg, barcc, IdivRg, numframes, Imin, masscut, field);
% toc

%% end