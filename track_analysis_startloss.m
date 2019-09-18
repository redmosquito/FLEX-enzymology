%% track_analysis_startloss
% This code will calculate take the histogram of trajectory lengths and fit
% the portion from the minimum time for flow extension through the stop of
% flow for a cell with 20 min flow-on time (2800 frames (23.3 mins) (2400 frames of flow after 400
% frames of no flow0
% in order to estimate quantitatively when the bead loss rate hits its
% maximum
%
% Need to be in parent directory for trajectories. Inputs are fovn and ver
%
% Drew Drabek
% Blacklow and Loparo Labs
% Harvard Medical School
% July 2018

%% load the data needed for determining the start time

% open the "TrackAnalysis" file and the lengths file
%   we need to have the extend and extlen variables and the length variable
    
    % example for AS11-169 flow extension
        % TrackAnalysis: load('\\molecule\blacklowlab\HomeDirectories\adrabek\Microscopy\BABarcus\AD10\AD10-169\chamber1\adig-s-Pmag9mm-cutting-A17cat-1uM_1\Track_sorting_fov169_v1\TrackAnalysis_fov_169_v1.mat')
        % TrackAnalysis_params: load('\\molecule\blacklowlab\HomeDirectories\adrabek\Microscopy\BABarcus\AD10\AD10-169\chamber1\adig-s-Pmag9mm-cutting-A17cat-1uM_1\Track_sorting_fov169_v1\TrackAnalysis_params_fov_169_v1.mat')
        % lengths: load('\\molecule\blacklowlab\HomeDirectories\adrabek\Microscopy\BABarcus\AD10\AD10-169\chamber1\adig-s-Pmag9mm-cutting-A17cat-1uM_1\Bead_tracking\lengths_fov169.mat')

%% find the lengths of interest
% we want to fit only the lengths that occur between the end of the
% flow-extension procedure and the stop of protease flow

% sort the length file
clear short; clear long; clear just;

    short = find(lengths <= 2800);     % the longest acceptable
    long = find(lengths >= (extend +extlen)); % the shortest acceptable
    just = intersect(short,long);   % the intersection of short and long
    
% fit the tracks from the just indices to a normal distribution
    fit = fitdist(lengths(just),'normal'); 
            % define the x and y data for the plot
            XLim = [(extend+extlen) 2800];
            XLim = XLim + [-1 1] * 0.01 * diff(XLim);
            XGrid = linspace(XLim(1),XLim(2),100);
            YPlot = pdf(fit,XGrid);
            
    % extract the mean and standard deviateion (sd)
    mu=fit.mu
    sd=fit.sigma 
    
    % assign the start of loss to one standard deviation before the mean
    startloss=mu-sd

%% plot the result
  
% overlay the histogram with the normal distribution
    figure
    hold on
    
    histogram(lengths(just),'BinWidth',100,'Normalization', 'pdf')
    
    plot(XGrid,YPlot,'LineWidth',2,'Color',[0.65 0.65 0.65]);
 
           % include legend with number of particles and mean/SEM fit
            legend({[' n = ' num2str(length(lengths))];...
                [' Normal fit (\mu = ' num2str(mu,4) ...
                ', \sigma = ' num2str(sd,3) ')']},...
            'Location','northeast')
            set(gca,'FontSize',10)
            set(gca,'box','off')
            xlabel('frame number')
            ylabel('probability density')
        hold off
    

% end