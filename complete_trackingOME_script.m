%% Complete tracking script
% Feature Finding mpretrack_init Template
% Tracking fancytrack Template
% MSD calculation MSD Template
% 
% Based on particle tracking code from http://people.umass.edu/kilfoil/tools.php
% Maria Kilfoil at UMass Amherst
%
% Drew Drabek
% Blacklow and Loparo Labs, Harvard Medical School, BCMP

function complete_trackingOME_script(fovn)
tic
%% mpretrack_init_template 
% Set Parameters for pretracking
basepath = [cd,'\'];    % Document the working directory for the images
   % basepath = [cd,'/'];    % Document the working directory for mac
    
%% a few constants that can change based on experimental setup
          %set the frame rate, number frames per second
        fr=2; % two is standard, can be higher for fast movies
%         fr=1;   
%         fr=4;
%          fr=6;
        
        %document the binning
        bin=1; % 1 is standard, can be 2, 4, 8
        
        % set the conversion from pixels to microns
%         conv=1.6125*bin;  %for 4x
        conv=0.6450*bin;  %for 10x
%         conv=0.16125*bin; %for 40x
        
        
%% Run the script mpretrack_renameFOV.m
%   the goal is to rename image files with fov prefix that matches the
%   folder name they reside in. Script works on tiff files saves renamed
%   and converted images in a subdirectory

renameFOVBF(fovn,'tif');

%% Define your parameters
% permissive 7,0,100,1,0
% strict 5,2000,40,0.2,0

% parameters for 10x tracking
featsize = 5;   %radius of reatures in pixels, The size of the feature you want to find.
barint = 2000;   %The minimum intensity you want to accept.
barrg = 45;      %The maximum Rg squared you want to accept.
barcc = 0.2;    %The maximum eccentricity you want to accept.
IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
%fovn = 1;    %ID# for the series of images (typically one field of view
frame = 0;      %ID# for the frame in the series
Imin=50;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
masscut=0;      %(optional) the masscut parameter for feature2D to remove false positives before rifining the position to speed up the code.
field=2;        %(optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
                %    image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
                %    set to 2 for progressive scan cameras. Defaults to 2.
                
% % parameters for 40x tracking (with binning)
% featsize =10;   %radius of reatures in pixels, The size of the feature you want to find.
% barint = 500;   %The minimum intensity you want to accept.
% barrg = 30;      %The maximum Rg squared you want to accept.
% barcc = 0.7;    %The maximum eccentricity you want to accept.
% IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
% %fovn = 1;    %ID# for the series of images (typically one field of view
% frame = 0;      %ID# for the frame in the series
% Imin=50;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
% masscut=0;      %(optional) the masscut parameter for feature2D to remove false positives before rifining the position to speed up the code.
% field=2;        %(optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
%                 %    image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%                 %    set to 2 for progressive scan cameras. Defaults to 2.

% % parameters for 60x tracking of fluorescent spots
% featsize =4;   %radius of reatures in pixels, The size of the feature you want to find.
% barint = 40000;   %The minimum intensity you want to accept.
% barrg = 12;      %The maximum Rg squared you want to accept.
% barcc = 0.5;    %The maximum eccentricity you want to accept.
% IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
% %fovn = 1;    %ID# for the series of images (typically one field of view
% frame = 0;      %ID# for the frame in the series
% Imin=50;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
% masscut=0;      %(optional) the masscut parameter for feature2D to remove false positives before rifining the position to speed up the code.
% field=2;        %(optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
%                 %    image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%                 %    set to 2 for progressive scan cameras. Defaults to 2.

% % parameters for binned 10x tracking of 3um beads
% featsize =4;   %radius of reatures in pixels, The size of the feature you want to find.
% barint = 400;   %The minimum intensity you want to accept.
% barrg = 20;      %The maximum Rg squared you want to accept.
% barcc = 0.15;    %The maximum eccentricity you want to accept.
% IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
% %fovn = 1;    %ID# for the series of images (typically one field of view
% frame = 0;      %ID# for the frame in the series
% Imin=50;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
% masscut=0;      %(optional) the masscut parameter for feature2D to remove false positives before rifining the position to speed up the code.
% field=2;        %(optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
%                 %    image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%                 %    set to 2 for progressive scan cameras. Defaults to 2.


% % parameters for 4x tracking
% featsize = 2;   %radius of reatures in pixels, The size of the feature you want to find.
% barint = 1000;   %The minimum intensity you want to accept.
% barrg = 4.5;      %The maximum Rg squared you want to accept.
% barcc = 0.15;    %The maximum eccentricity you want to accept.
% IdivRg = 0;     %Minimum ratio of Intensity/pixel to be accepted (integratedintensity / Rg squared of feature)
% %fovn = 1;    %ID# for the series of images (typically one field of view
% frame = 0;      %ID# for the frame in the series
% Imin=50;         %(optional) the minimum intensity for a pixel to be considered as a potential feature.
% masscut=0;      %(optional) the masscut parameter for feature2D to remove false positives before rifining the position to speed up the code.
% field=2;        %(optional) set to 0 or 1 if image is actually odd or even field of an interlaced 
%                 %    image. All the masks will then be constructed with a 2:1 aspect ratio. Otherwise 
%                 %    set to 2 for progressive scan cameras. Defaults to 2.

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

% save the mpretrack_init figure as a png file
openfig(fullfile(basepath, 'mpretrack_init/', ['fov' num2str(fovn) '_pretrackfig.fig']));
    print('-f',fullfile(basepath,'mpretrack_init',['fov' num2str(fovn) '_pretrackfig.png']),'-dpng');
close(gcf);

%% Feature finding loop function for all fields of view
%calculate the number of frames
%numframes=length(dir(fullfile(['fov',num2str(fovn)],'*.tif')))-1;
load(fullfile(['fov' num2str(fovn) '_times.mat'])) 
numframes=length(time)-1;
%numframes=600; %uncomment this line to test tracking for fewer frames

%pretrack all the frames specified in numframes
tic
mpretrack(basepath, fovn, featsize, barint, barrg, barcc, IdivRg, numframes, Imin, masscut, field);
toc

%% fancytrack_MSD_template
% Drew Drabek Blacklow Lab Harvard Medical School BCMP
% 20161129

% run tracking algorithm with pretrack as input
% Workhorse for tweezers 15, 30, 1  ( works for most cases)
% Workhorse for tweezers 18, 300, 5 (alternative for fast extension)
% Permissive 10,10,10
% Strict 1,100,2
% Fluorescent 2, 2, 2
maxdisp=22;     % number of steps that can be taken between frames
goodenough=99;  % minimum length of track
memory=4;      % number of frames accepted without particle present

%run fancytrack
tic
fancytrack(basepath,fovn,featsize,maxdisp,goodenough,memory)
toc

%save the parameters used for fancytrack
tic
save(fullfile(basepath,'Bead_tracking/res_files',['fov' num2str(fovn) '_pretrackparams']));
load(fullfile(basepath, 'Bead_tracking/res_files', ['res_fov' num2str(fovn) '.mat']));
csvwrite(fullfile(basepath, 'Bead_tracking/res_files', ['res_fov' num2str(fovn) '.csv']),res); 
toc

%% analysis of tracking and saving variables
% M2/MT matrix format:
% 1 row per bead per frame, sorted by bead ID then by frame number.
% columns are:
% 1:2 - X and Y positions (in pixels) xpos ypos
% 3   - Integrated intensity int
% 4   - Rg squared of feature rg
% 5   - eccentricity cc
tic
%M2_header xpos ypos int rg cc 
csvwrite(fullfile(basepath, 'mpretrack_init', ['M2_fov' num2str(fovn) '.csv']),M2);

%MT_header xpos ypos int rg cc 
csvwrite(fullfile(basepath, 'Feature_finding', ['MT_fov' num2str(fovn) '.csv']),MT);

% res matrix format:
% 1 row per bead per frame, sorted by bead ID then by frame number.
% columns are:
% 1:2 - X and Y positions (in pixels) xpos ypos
% 3   - Integrated intensity int
% 4   - Rg squared of feature rg
% 5   - eccentricity cc
% 6   - frame # 
% 7   - time of frame
% 8   - Bead ID

%res_header xpos ypos int rg cc num time beadID
load(fullfile(basepath, 'Bead_tracking/res_files', ['res_fov' num2str(fovn) '.mat']));
csvwrite(fullfile(basepath, 'Bead_tracking/res_files', ['res_fov' num2str(fovn) '.csv']),res); 

toc
%%  Interface with @msdanalyzer
%   This section writes the output of the particle tracking for each bead
%   into a cell array called "tracks" so that the individual MSD for each
%   particle can be calculated and plotted
%
%   REFERENCE 
%   "Mean square displacement analysis of particle trajectories."
%
%   https://tinevez.github.io/msdanalyzer/
%
%   Nadine Tarantino, Jean-Yves Tinevez, Elizabeth Faris Crowell, 
%   Bertrand Boisson, Ricardo Henriques, Musa Mhlanga, Fabrice Agou, 
%   Alain Israël, and Emmanuel Laplantine. TNF and IL-1 exhibit distinct 
%   ubiquitin requirements for inducing NEMO-IKK supramolecular structures.
%   J Cell Biol (2014) vol. 204 (2) pp. 231-45
%
%   How to use 
%       in interface with Kilfoil's Matlab adaptation of Cocker and Weeks
%       code written in IDL

%   Need to have opened a file with linked particle trajectories. In this
%   case this is the output of the function fancytrack.m which outputs an
%   n x 8 double with n particle locations total from an image sequence.
%   The columns 1 and 2 are the x and y positions, respectively
%   The column 7 is the timestamp for the the position
%   The column 8 is the particle ID
%
%   This code will use this res file and for each unique particle ID, will
%   write the timestamp, x-position, and y-position to a n x 3 double in a
%   cell in the cell array tracks
%
%   This array tracks is in the format to be recognized by the @msdanalyzer
%   code for ease of analysis of particles for displacement and
%   identification of particles that meet criteria to be identified as
%   specific tethers

% Transfer tracks to msdanalyzer, plot tracks and msd analysis

    % write the empty tracks cell arrays to be generated
    tracks={};
    tracksum={};
    tracksres={};
    tracksresum={};

     %conv=1.6125;  %for 4x
     %conv=0.6450;  %for 10x
     %conv=0.16125; %for 40x
    
    % write the empty frames cell array and double to be generated
        frames={};
        framesct=[];
    
    tic
    % use a loop to iterate through the number of particles and save the
    % output in bot pixels (tracks and tracksres) and microns (+um suffix)
    for i = 1:max(res(:,8))
        % assign the indices for the ith particle
        indices=find(res(:,8)==i);

        % for msdanalyzer write the time(col7) xpos(col1) and ypos(col2) to the ith tracks cell
        tracks{i}=[res(indices(1):indices(end),7) res(indices(1):indices(end),1) res(indices(1):indices(end),2)];
        tracksum{i}=[res(indices(1):indices(end),7) res(indices(1):indices(end),1).*conv res(indices(1):indices(end),2).*conv];
        
        % for saving res matrix as a cell array
        tracksres{i}=[res(indices(1):indices(end),1) res(indices(1):indices(end),2) res(indices(1):indices(end),3) res(indices(1):indices(end),4).*(conv^2) res(indices(1):indices(end),5) res(indices(1):indices(end),6) res(indices(1):indices(end),7) res(indices(1):indices(end),8)];
        tracksresum{i}=[res(indices(1):indices(end),1).*conv res(indices(1):indices(end),2).*conv res(indices(1):indices(end),3) res(indices(1):indices(end),4) res(indices(1):indices(end),5) res(indices(1):indices(end),6) res(indices(1):indices(end),7) res(indices(1):indices(end),8)];
        
    end
    
    clear i indices
    load(fullfile(basepath, 'Feature_finding/', ['MT_' num2str(fovn) '_Feat_Size_' num2str(featsize) '.mat']));
    
    % use a loop to save the number of accepted objects per frame
    for i = 1:max(MT(:,6))
        % assign the indices for the ith particle
        indices=find(MT(:,6)==i);
        if length(indices) > 1
        
        % for MT file save the number of objects per frame
        frames{i}=[MT(indices(1):indices(end),7) MT(indices(1):indices(end),1) MT(indices(1):indices(end),2)];
        
        framesct(i,1)=MT(indices(1),6); %frame number
        framesct(i,2)=length(indices);  %number of particles 
        else
            continue
        end
        
        
    end
    
    tracks=tracks';
    tracksres=tracksres';
    tracksum=tracksum';
    tracksresum=tracksresum';

    % save tracks file
    save(fullfile(basepath,'Bead_tracking',['tracks_fov' num2str(fovn)]),'tracks' , 'tracksum');
    
    % save tracksres file
    save(fullfile(basepath,'Bead_tracking',['tracksres_fov' num2str(fovn)]),'tracksres' , 'tracksresum');
    toc
    
    tic
    % Initialize ma analyzer
    SPACE_UNITS = 'pixels';
    TIME_UNITS = 's';
    ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

    % Load tracks into ma class
    ma = ma.addAll(tracks);
    
    % output of ma to command window
    disp(ma)
  

    % save a figure with the tracks
    figure;
            T=['fov' num2str(fovn) ' Particle Tracks'];
            ma.plotTracks;
            ma.labelPlotTracks;
            set(gca,'ydir','reverse')% invert the y-axis to make consistent with mpretrack
            savefig(fullfile(basepath, 'Bead_tracking',['PlotTracks_' SPACE_UNITS '_fov' num2str(fovn) '.fig']));
            print('-f',fullfile(basepath, 'Bead_tracking',['PlotTracks_' SPACE_UNITS '_fov' num2str(fovn) '.png']),'-dpng');
    close(gcf);
    
    % Initialize ma analyzer
    SPACE_UNITS = 'microns';
    TIME_UNITS = 's';
    ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

    % Load tracks into ma class
    ma = ma.addAll(tracksum);
    
    % output of ma to command window
    disp(ma)
  

    % save a figure with the tracks
    figure;
            T=['fov' num2str(fovn) ' Particle Tracks'];
            ma.plotTracks;
            ma.labelPlotTracks;
            set(gca,'ydir','reverse')% invert the y-axis to make consistent with mpretrack
            savefig(fullfile(basepath, 'Bead_tracking',['PlotTracks_' SPACE_UNITS '_fov' num2str(fovn) '.fig']));
            print('-f',fullfile(basepath, 'Bead_tracking',['PlotTracks_' SPACE_UNITS '_fov' num2str(fovn) '.png']),'-dpng');
    close(gcf);
 
% plot crude analysis (akin to 'ImageJ' analysis, but eliminated clusters)
    figure;
        hold on
        T=['fov' num2str(fovn) ' Crude # Particles vs Time'];
            xlabel('time, s');
            ylabel('fraction of initial objects remaining');
            ymax=max(framesct(:,2)); %maximum for count
            ylim([0 1]);
            plot(framesct(:,1),framesct(:,2)./ymax);
        hold off 
            savefig(fullfile(basepath, 'Bead_tracking',['PlotFrac_raw_fov' num2str(fovn) '.fig']));
            print('-f',fullfile(basepath, 'Bead_tracking',['PlotFracvs_raw_fov' num2str(fovn) '.png']),'-dpng');
    close(gcf);
    clear fr;
    toc
%% Track length distribution
%
    %calculate a histogram of the frequency of different track lengths
    
    % output length of each trajectory from tracks into lengths
    lengths=cellfun(@length,tracksres);
    save(fullfile(basepath,'Bead_tracking',['lengths_fov' num2str(fovn)]),'lengths');
    
    %plot the histogram 
    figure;
        hist(lengths,round(sqrt(length(lengths))));            
        T=['fov' num2str(fovn) ' Track Length Histogram'];
            xlabel('track length, s');
            ylabel('number of tracks');
        savefig(fullfile(basepath, 'Bead_tracking',['LengthHistogram_fov' num2str(fovn) '.fig']));
        print('-f',fullfile(basepath, 'Bead_tracking',['LengthHistogram_fov' num2str(fovn) '.png']),'-dpng');

    close(gcf);
    
%% count the number of beads bound taken from still images


    
%% APPENDIX - Code not used in script version of this file

% calculate mean square displacement and analyze it
% %open res tarcking file
% load(fullfile('Bead_tracking','res_files',['res_fov' num2str(fovn) '.mat'])) 
% 
% %create minimal matrix to calculate MSD
% linked=res(:,[1:2,7:8]);    
% 
% %calculate MSD in x and y
% tic
% msd=MSD(linked);
% toc
% % the MSD will have three columns:
% %   Column 1: lag time (in frames)
% %   Column 2: MSD (in pixels)
% %   Column 3: number of observations in average
% 
% %MSD_header lagtime MSD numobs
% csvwrite(fullfile(basepath, 'Bead_tracking', ['msd_fov' num2str(fovn) '.csv']),msd);
%    
% %save plot of MSD vs lag time
% figure;
%     loglog(msd(:,1),msd(:,2));
%         T=['fov' num2str(fovn) ' Mean Square Displacement'];
%         xlabel('lag time in frames');
%         ylabel('msd in pixels');
%     savefig(fullfile(basepath, 'Bead_tracking',['MSD_fov' num2str(fovn) '.fig']));
%     print('-f',fullfile(basepath, 'Bead_tracking',['MSD_fov' num2str(fovn) '.png']),'-dpng');
% 
% close(gcf);

% % % save scatter plot of tracked particles
% % figure;    
% % % scatter(res(:,1),res(:,2),res(:,4),res(:,8),'fill'); % size prop Rg
% % scatter(res(:,1),res(:,2),1,res(:,8),'fill'); % size = 1
% %     T=['fov' num2str(fovn) ' Tracking Scatter Plot'];
% %     title(T);
% %     xlabel('X (pixels)');
% %     ylabel('Y (pixels)');
% %     savefig(fullfile(basepath, 'Bead_tracking',['Bead_tracking_fov' num2str(fovn) '.fig']));
% %     print('-f',fullfile(basepath, 'Bead_tracking',['BeadTracking_fov' num2str(fovn) '.png']),'-dpng');
% % 
% % close(gcf);
% 
% % % %% Interface with MSD analyzer - Calculate MSD and mean MSD
% % % %calculate and save a figure of the msd
% % %     tic
% % %     ma = ma.computeMSD;
% % %     toc
% % %     ma.msd;
% % % 
% % %     %plot and save MSD tracks, MSD vs delay time
% % %     figure;
% % %             T=['fov' num2str(fovn) ' Mean Square Displacement'];
% % %             ma.plotMSD;
% % %             savefig(fullfile(basepath, 'Bead_tracking',['MSDanalyzer_fov' num2str(fovn) '.fig']));
% % %             print('-f',fullfile(basepath, 'Bead_tracking',['MSDanalyzer_fov' num2str(fovn) '.png']),'-dpng');
% % % 
% % %     close(gcf);
% % % 
% % %     %Plot mean MSD curves
% % %     figure;
% % %             T=['fov' num2str(fovn) ' Mean Mean Square Displacement'];
% % %         cla
% % %         ma.plotMeanMSD(gca, true)        
% % %         savefig(fullfile(basepath, 'Bead_tracking',['MeanMSDanalyzer_fov' num2str(fovn) '.fig']));
% % %         print('-f',fullfile(basepath, 'Bead_tracking',['MeanMSDanalyzer_fov' num2str(fovn) '.png']),'-dpng');
% % % 
% % %     close(gcf);
% % % 
% % % % Track length distribution
% % %     %calculate a histogram of the frequency of different track lengths
% % %     
% % %     % output length of each trajectory from tracks into lengths
% % %     lengths=cellfun(@length,tracks);
% % %     save(fullfile(basepath,'Bead_tracking',['lengths_fov' num2str(fovn)]),'lengths');
% % %     %plot the histogram
% % %      
% % %     figure;
% % %         hist(lengths,50);            
% % %         T=['fov' num2str(fovn) ' Track Length Histogram'];
% % %             xlabel('track length, s');
% % %             ylabel('number of tracks');
% % %         savefig(fullfile(basepath, 'Bead_tracking',['LengthHistogram_fov' num2str(fovn) '.fig']));
% % %         print('-f',fullfile(basepath, 'Bead_tracking',['LengthHistogram_fov' num2str(fovn) '.png']),'-dpng');
% % % 
% % %     close(gcf);
% % %     
% % % % MSD Histogram
% % % 
% % % % define double to contain msd values
% % %     msdhist=[];
% % %     dispxy={};
% % %     % calculate the displacements for each particle
% % %     for j=1:length(ma.tracks)
% % %         
% % %         % iterate for each particle for the length of the track
% % %         for i=2:length(ma.tracks{j,1});
% % %         
% % %         %displacement (x-xi) + (y-yi)
% % %         dispxy{j,1}(i,1)=((ma.tracks{j,1}(i-1,2)-ma.tracks{j,1}(i,2))+(ma.tracks{j,1}(i-1,3)-ma.tracks{j,1}(i,3)));
% % %         end
% % %     end
% % %     
% % %     save(fullfile(basepath,'Bead_tracking',['displacements_fov' num2str(fovn)]),'dispxy');
% % % 
% % %     %plot a histogram of all the displacements
% % %     figure;
% % %         hist(cell2mat(dispxy(:,1)),100)
% % %         T=['fov' num2str(fovn) ' Displacements Histogram '];
% % %             xlabel('displacement, pixels');
% % %             ylabel('number of measurments');
% % %         savefig(fullfile(basepath, 'Bead_tracking',['DisplacementHistogram_fov' num2str(fovn) '.fig']));
% % %         print('-f',fullfile(basepath, 'Bead_tracking',['DisplacementHistogram_fov' num2str(fovn) '.png']),'-dpng');
% % % 
% % %     close(gcf);
% % %         
% 
% % % % calculate MSD for each particle and populate msdhist
% % %     for i=1:length(ma.msd)
% % %         msdx=ma.msd{i,1};
% % %         msdy=ma.msd{i,2};
% % %         msdhist(i,1)=fitlm(ma.msd{i,1}(:,1),ma.msd{i,1}(:,2));
% % %     end
% % %     save(fullfile(basepath,'Bead_tracking',['indivMSDx_fov' num2str(fovn)]),'msdx');
% % %     save(fullfile(basepath,'Bead_tracking',['indivMSDy_fov' num2str(fovn)]),'msdy');
% 
% 
% % % % Particle Velocities
% % %     
% % %     % displacements computed over varying time interval
% % %     % The velocities are returned in the shape of one 3D double array, arranged as for MSD curves: [T Vx Vy ...] 
% % %     v = ma.getVelocities %#ok<NOPTS>
% % %     
% % %         %Let's pool all particles together:
% % % 
% % %     V = vertcat( v{:} );

    
toc
%% end
% Drew Drabek   ver1   20161129
% Drew Drabek   ver2   20170605
% Drew Drabek   ver3   20170613