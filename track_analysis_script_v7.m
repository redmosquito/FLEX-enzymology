%% Track Analysis - Detailed analysis of particle motion
% This script can be run after a set of trajectories has been created for a
% movie of tethered colloids is taken. First particles are tracked with the
% complete_tracking script and then they are loaded into the msdanalyzer
% matalb code and the "tracks" cell array containing time stamp, x and y
% positions is the input for this analysis
%   
%   fovn    = field of view number
%   tracks  = describe the set of trajectories to be analyzed
%               i.e. AllTracks or HighlyMobileTracks
%   id      = unique identifier for this analysis
%   start   = particle ID of first track to subplot in singletracks
%   stop    = particle ID of last track to subplot in singletracks
%   initial = frame to start averaging for initial position calculation
%   extend  = frame to start averaging for extended position calculation
%   initlen = number of frame to average for displacement initial position
%   extlen  = number of frame to average for displacement extended position
%   toolegit= lower limit for legit trajectories
%           = upper limit for mediocre tracks
%   toquit  = upper limit for rejected tracks
%           = lower limit for mediocre
%   freakish= upper limit for legit trajectories
%   conv    = the conversion microns per pixel
% 
%       % recommended values for thresholds with 10x objective and 8um tether
%           %toolegit is 7.5 to toquit 2 for 1pN
%               freakish=13.5;
%               toolegit=7.5;
%               toquit=2;
%           %toolegit 3 to toquit 1 for 3.6pN
%               freakish=11;
%               toolegit=6;
%               toquit=1.5;
%           %toolegit 2 to toquit 0.5 for 5.4pN.
%               freakish=8;
%               toolegit=3.5;
%               toquit=1;
%
% recommended values for thresholds with 4x objective and 8um tether
%           %toolegit is 7.5 to toquit 2 for 1pN
%               freakish=4.5;
%               toolegit=3.25;
%               toquit=1;
%
%    conv=0.16125; %for 40x
%    conv=0.6450; % for 10x
%    conv=1.6125; % for 4x

%   Drew Drabek
%   Blacklow and Loparo Labs
%   BCMP, Harvard Medical School

function track_analysis_script_v7(fovn,id,ver,start,initial,initlen,...
    extend,extlen,toquit,toolegit,freakish,conv)
%% Generate a sampling of trajectories focusing on x and y positions

% fovn=491;

% redefine basepath to load and save variables
    %basepath = [cd,''];    % Document the working directory for the images
     basepath = [cd,'/'];    % Document the working directory for the images


% Load tracking parameters for use of variables
        load(fullfile([basepath '/Bead_tracking/res_files/fov' num2str(fovn) '_pretrackparams.mat']))
        basepath = [cd,'/'];    % Document the working directory for the images
        load(fullfile([basepath '/Bead_tracking/res_files/fov' num2str(fovn) '_pretrackparams.mat']),'maxdisp','memory')

    %load(fullfile([basepath 'Bead_tracking/res_files/fov' num2str(fovn) '_pretrackparams.mat']))

    % basepath = [cd,''];    % Document the working directory for the images
     basepath = [cd,'/'];    % Document the working directory for the images


% define the trajectories of interest
    load(fullfile(basepath, '/Bead_tracking/', ['tracksres_fov' num2str(fovn) '.mat']));
        %load(fullfile(basepath, 'Bead_tracking/', ['tracksres_fov' num2str(fovn) '.mat']));

basepath = [cd,'/'];    % Document the working directory for the images

% save the input parameters
    mkdir(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)]);
    save(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
        ['TrackAnalysis_params_fov_' num2str(fovn) '_v' num2str(ver) '.mat']),...
        'fovn','id','ver','initial','initlen','extend','extlen','toquit',...
        'toolegit','freakish','basepath');
    

%% initialize for mac
% 
% % redefine basepath to load and save variables
%     basepath = [cd,'/'];    % Document the working directory for the images
% 
% 
% % Load tracking parameters for use of variables
%         load(fullfile([basepath 'Bead_tracking/res_files/fov' num2str(fovn) '_pretrackparams.mat']))
%             
%         basepath = [cd,'/'];    % Document the working directory for the images
% 
%         load(fullfile([basepath 'Bead_tracking/res_files/fov' num2str(fovn) '_pretrackparams.mat']),'maxdisp','memory')
% 
%     %load(fullfile([basepath 'Bead_tracking/res_files/fov' num2str(fovn) '_pretrackparams.mat']))
% 
%      basepath = [cd,'/'];    % Document the working directory for the images
% 
% 
% % define the trajectories of interest
%     load(fullfile(basepath, 'Bead_tracking/', ['tracksres_fov' num2str(fovn) '.mat']));
%         %load(fullfile(basepath, 'Bead_tracking/', ['tracksres_fov' num2str(fovn) '.mat']));
% 
% basepath = [cd,'/'];    % Document the working directory for the images
% 
% % save the input parameters
%     save(fullfile(basepath,'Bead_tracking/',...
%         ['TrackAnalysis_params_fov' num2str(fovn) '.mat']),...
%         'fovn','id','initial','extend','initlen','extlen');
%     
%% load the tracks to analyze and save a name
% save the name of the tracks to be analyzed into a character array
    name=whos('tracksres');
    name=name.name;

    % ensure that tracksres array is ordered correctly and save it into
    % trackplot variable (tracksres can be substituted for another array if
    % desired)
    sz=size(tracksres);
    if sz(1,1) < sz(1,2)
        tracksres=tracksres';
        trackplot=tracksres;
        'if'
    else
        trackplot=tracksres;
        fprintf('else');
    end
    
%     % manual setting parameters in the function
%         %fovn=1630;
%         id='60sInPmag13mm';
%         ver=2;
%         start=151;
%         initial=1079;
%         initlen=120;
%         extend=1301;
%         extlen=120;
%         toquit= 2;
%         toolegit=7.5;
%         freakish=13.5;
    
    %standard analysis is to plot 25 tracks
    stop=start+24;
    
    % concatenate the name identifier for the saved variables and files
    name = [id '_' num2str(maxdisp) 'step_' num2str(memory) 'mem_' num2str(initial) '-' num2str(extend) '_bins' ];
 
    
%   name=whos('tracksreslegitgated');
%   name=name.name;
%   trackplot=tracksreslegitgated;
%
%         start=301;
%         stop=325;
    
%% plot the x-position data
    % plot x position vs time normalized
    name1=['x-position vs time ' num2str(start) ' to ' num2str(stop)];
    figure('Name',[name1 ' norm']);
        
        i=0;
        for j=start:stop
            i=i+1;
            subplot(5,5,i)
            % %plot the position vs time in the ith subplot
            % plot(trackplot{j,1}(1:end-10,7),trackplot{j,1}(1:end-10,1),'r--o');
            
            %normalize positions  ((divide by minimum subtract 1 )x100)
            normx=trackplot{j,1}(1:end-10,1)-min(trackplot{j,1}(1:end-10,1));
            normy=trackplot{j,1}(1:end-10,2)-min(trackplot{j,1}(1:end-10,2));
            time=trackplot{j,1}(1:end-10,7);
            
            %plot the position vs time in the ith subplot
            plot(time,normx,'r-');
%             hold on
% %                 f1=fit(time,normx,'smoothingspline');
% %                     plot(f1,'k--')
%                 f2=fit(time,normx,'linearinterp');
%                     plot(f2,'r--')
%             hold off
            
            subplot(5,5,i)
                %ymin=1;
                %ymax=5;
                %ylim([ymin ymax]);
                %xlim([0 400])
                title(['particle' num2str(j)]);

        end
        clear i;
        
    save(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
        ['SingleTracks_XposVsTime_' num2str(start) '-' num2str(stop) name '_norm.fig']));
    print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
        ['SingleTracks_XposVsTime_' num2str(start) '-' num2str(stop) name '_norm']),'-dpng');
   
    close(gcf);
    
    % Plot the x-position data raw (no normalization)
        figure('Name',[name1 ' raw']);
        
        i=0;
        for j=start:stop
            i=i+1;
            subplot(5,5,i)
            hold on
            %plot the position vs time in the ith subplot
            scatter(trackplot{j,1}(1:end-10,7),trackplot{j,1}(1:end-10,1),50,trackplot{j,1}(1:end-10,7),'.');
            line(trackplot{j,1}(1:end-10,7),trackplot{j,1}(1:end-10,1),'Color','r')

            % % plot the position vs time in the ith subplot
            %  plot(time,normx,'r-',time,normy,'k-');
            
            subplot(5,5,i)
%                 ymin=mean(trackplot{j,1}(1:end-10,1))-10;
%                 ymax=mean(trackplot{j,1}(1:end-10,1))+10;
%                 ylim([ymin ymax]);
                %xlim([0 400])
                title(['particle' num2str(j)]);
            hold off
    
        end
        clear i;
        
        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['SingleTracks_XposVsTime_' num2str(start) '-' num2str(stop) name '_raw.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['SingleTracks_XposVsTime_' num2str(start) '-' num2str(stop) name '_raw']),'-dpng');
   
    close(gcf);
    
%% plot the y-position data
    % plot y position vs time normalized
    name2=['y-position vs time ' num2str(start) ' to ' num2str(stop)];
    figure('Name',[name2 ' norm']);
        
        i=0;
        for j=start:stop
            i=i+1;
            subplot(5,5,i)
            
            %normalize positions  ((divide by minimum subtract 1 )x100)
            normy=trackplot{j,1}(1:end-10,2)-min(trackplot{j,1}(1:end-10,2));
            time=trackplot{j,1}(1:end-10,7);
            
            %plot the position vs time in the ith subplot
            plot(time,normy,'b-');
            
            subplot(5,5,i)
                %ymin=1;
                %ymax=5;
                %ylim([ymin ymax]);
                %xlim([0 400])
                title(['particle' num2str(j)]);

        end
        clear i;
        
        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['SingleTracks_YposVsTime_' num2str(start) '-' num2str(stop) name '_norm.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['SingleTracks_YposVsTime_' num2str(start) '-' num2str(stop) name '_norm']),'-dpng');
   
    close(gcf);
    
        % Plot the raw y-position data (no normalization)
            figure('Name',[name2 ' raw']);

            i=0;
            for j=start:stop
                i=i+1;
                subplot(5,5,i)
            hold on
            %plot the position vs time in the ith subplot
            scatter(trackplot{j,1}(1:end-10,7),trackplot{j,1}(1:end-10,2),50,trackplot{j,1}(1:end-10,7),'.');
            line(trackplot{j,1}(1:end-10,7),trackplot{j,1}(1:end-10,2),'Color','b')

                    %ymin=mean(trackplot{j,1}(:,1))-10;
                    %ymax=mean(trackplot{j,1}(:,1))+10;
                    %ylim([ymin ymax]);
                    %xlim([0 400])
                    title(['particle' num2str(j)]);
            hold off
            end
            clear i;

            savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
                ['SingleTracks_YposVsTime_' num2str(start) '-' num2str(stop) name '_raw.fig']));
            print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
                ['SingleTracks_YposVsTime_' num2str(start) '-' num2str(stop) name '_raw']),'-dpng');

        close(gcf);
        
%% plot x- and y-position vs time normalized
    
    % plot the normalized y-position data
        name3=['x- and y-position vs time ' num2str(start) ' to ' num2str(stop)];
        figure('Name',[name3 ' norm']);

            i=0;
            for j=start:stop
                i=i+1;
                subplot(5,5,i)
               

                %normalize positions  ((divide by minimum subtract 1 )x100)
                normx=trackplot{j,1}(1:end-10,1)-min(trackplot{j,1}(1:end-10,1));
                normy=trackplot{j,1}(1:end-10,2)-min(trackplot{j,1}(1:end-10,2));
                time=trackplot{j,1}(1:end-10,7);

                %plot the position vs time in the ith subplot
                plot(time,normx,'r-',time,normy,'b-');
%                 % plot x vs y normalized
%                 plot(normx,normy,20,'filled');

                subplot(5,5,i)
                    %ymin=1;
                    %ymax=5;
                    %ylim([ymin ymax]);
                    %xlim([0 400])
                    title(['particle' num2str(j)])

            end
            clear i;

            savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
                ['SingleTracks_X&YposVsTime_' num2str(start) '-' num2str(stop) name '_norm.fig']));
            print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
                ['SingleTracks_X&YposVsTime_' num2str(start) '-' num2str(stop) name '_norm']),'-dpng');

        close(gcf);
    
        % plot the x and y trajectories
            figure('Name',[name3 ' raw']);

            i=0;
            for j=start:stop
                i=i+1;
                subplot(5,5,i)

                %plot the position vs time in the ith subplot
                hold on
                %plot the position vs time in the ith subplot
                scatter(trackplot{j,1}(1:end-10,1),trackplot{j,1}(1:end-10,2),50,trackplot{j,1}(1:end-10,7),'.');
                %xlim([0 400])
                    title(['particle' num2str(j)])
                hold off

            end
            clear i;

            savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
                ['SingleTracks_X&YposVsTime_' num2str(start) '-' num2str(stop) name '_raw.fig']));
            print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
                ['SingleTracks_X&YposVsTime_' num2str(start) '-' num2str(stop) name '_raw']),'-dpng');

        close(gcf);
        
%% Exclusion of aberrant tracks that are less than a minimum time

 % this section is commented out because the initial tracking algorithm
 % limits track length to 30s
    % Load tracking parameters for use of variables
    load(fullfile(basepath, 'Bead_tracking', ['lengths_fov' num2str(fovn) '.mat']));
    
    %define cell array for tracks that are long enough 
    tracksreslong={};

    %find indices for long tracks
    long=find(lengths>=(extend+extlen));
            
    % write legitimate elements of tracksres into the legit array
    tracksreslong=tracksres(long);
    
    %define cell array for particles starting at 0 (the gate)
    tracksresgated={};
    tracksreslonggated={};  
       
    %define cell array for particles that arrive after the gate
    tracksreslate={};
    tracksreslonglate={};   

    % iterate through the tracks and write particles starting at the gate
    % into tracksresgated cell array for all tracks
    for i = 1: length(tracksres)
       if tracksres{i}(1,7)==0
           tracksresgated{i,1}=tracksres{i};    %write gated tracks 
       else
           tracksreslate{i,1}=tracksres{i};     %write late starting tracks
       end
    end
    
    % iterate through the tracks and write particles starting at the gate
    % into tracksreslonggated cell array for all tracks
    for i = 1: length(tracksreslong)
       if tracksreslong{i}(1,7)==0
           tracksreslonggated{i,1}=tracksreslong{i};    %write gated tracks 
       else
           tracksreslonglate{i,1}=tracksreslong{i};     %write late starting tracks
       end
    end
    
    
    %% plot calculate x displacement

% tracks displacement matrix
trackcalc=tracksreslonggated;
tracksdisp=[length(trackcalc) 13];

% trackfit - cell array for the fitting of initial and extended
trackfit=cell([length(trackcalc) 8]);

%loop through each track, calculate the displacement
for i = 1:length(trackcalc)
     % determine if cell is empty, continue if true
    if isempty(trackcalc{i,1});
        continue
    else
        %calculate the minumum value of x and y for the ith track
        %exclude last 10 frames
        xmin=min(trackcalc{i}(1:end-10,1));
        ymin=min(trackcalc{i}(1:end-10,2));

        %calculate the maximum value of x and y for the ith track
        %exclude last 10 frames
        xmax=max(trackcalc{i}(1:end-10,1));
        ymax=max(trackcalc{i}(1:end-10,2));

        %calculate the max overall displacement for the ith track
        xdisp=xmax-xmin;
        ydisp=ymax-ymin;
        lindisp=sqrt((xmax-xmin)^2+(ymax-ymin)^2);

    %calculate the displacement due to flow with averaging from mean pos
        % average the extended position
        for j=1:extlen %iterate the extended position
            m=extend; %define the mth frame to start at for extended
            m=m+j;    %iterate to the jth frame after the start extended
            %calculate extended displacement as average of mth to m+jth frames
            xdispmean{i,1}(j,2)=trackcalc{i}(m,1);
            ydispmean{i,1}(j,2)=trackcalc{i}(m,2);
        end
        
        % average the initial position
        for k=1:initlen %iterate the initial position
            l=initial;  %define the lth frame to start at for initial
            l=l+k;      %iterate to the kth frame after the start, initial pos
            %calculate initial displacement as average of lth to l+kth frames
            xdispmean{i,1}(k,1)=trackcalc{i}(l,1);
            ydispmean{i,1}(k,1)=trackcalc{i}(l,2);
        end
        
    %calculate initial displacement from gaussian fitting of positions     
        % fit the scatterplot for the initial phase defined by initian and
        % meanlen (typically from 1-60 or 1-100frames, at 1s or 0.5s fr)
        fitin=fitgmdist(trackcalc{i}(initial:(initial+initlen),1:2),1);
            % assign initial position mean into first column
            din(1,1:2)=fitin.mu;        %mean x,y position
            din(2,1)=fitin.Sigma(1,1);   %variance xpos
            din(2,2)=fitin.Sigma(2,2);   %variance ypos
            
            % assign the initial fit into a cell array
            trackfit{i,1}=fitin;
            trackfit{i,3}=fitin.mu;    % means X,Y
            trackfit{i,4}=fitin.Sigma; % covariance matrix (XX, XY; YX, YY)
            trackfit{i,5}=sqrt(fitin.Sigma); % standard deviation
            
        %calculate extended displacement as average of lth to l+jth frames
        % fit the scatterplot for the extended phase defined by initian and
        % meanlen (typically from 200-250 or 400-500 frames at 1s or 0.5s fr)
        fitex=fitgmdist(trackcalc{i}(extend:extend+extlen,1:2),1);
            % assign extended position mean into first column
            dex(1,1:2)=fitex.mu;        %mean x,y position
            dex(2,1)=fitex.Sigma(1,1);   %variance xpos
            dex(2,2)=fitex.Sigma(2,2);   %variance ypos
            
            trackfit{i,2}=fitex;
            trackfit{i,6}=fitex.mu;    % means X,Y
            trackfit{i,7}=fitex.Sigma; % covariance matrix (XX, XY; YX, YY)
            trackfit{i,8}=sqrt(fitex.Sigma); % standard deviation
            
            dx=dex(1,1)-din(1,1);               % x disp
            sx= sqrt(dex(2,1) + din(2,1));      % x disp RMSE 
            
            dy=dex(1,2)-din(1,2);               % y disp
            sy= sqrt(dex(2,2) + din(2,2));      % y disp RMSE
            
            dl=sqrt((dex(1,1)-din(1,1))^2+(dex(1,2)-din(1,2))^2); %lin disp
            sl= sqrt(sy^2 + sx^2);       %lin disp RMSE, error propagated
            
        %     %uncomment this portion to do a spot check on the fitting
        %     figure;
        %         hold on
        %         scatter(trackcalc{i}(:,1),trackcalc{i}(:,2))
        %             %save the functions x1 and x2 into f
        %             fin=@(x1,x2)pdf(fitin,[x1 x2]);
        %             fex=@(x1,x2)pdf(fitex,[x1 x2]);
        %             % plot contours of f, output into f1
        %                 fc1=fcontour(fin,'c');
        %                 fc2=fcontour(fex,'m');
        %         hold off
        %         clear fin fex fc1 fc2
        %     close(gcf)

        % write the  displacements in x and y into the double tracksdisp
        %   column1 = maximum x displacement for entire track
        %   column2 = maximum y displacement for entire track
        %   column3 = maximum linear displacement for entire track
        %   column4 = x displacemnt between two time periods indicated
        %   column5 = y displacemnt between two time periods indicated
        %   column6 = linear displacement between two time periods indicated
        %   column7 = particle ID from column8 of trackcalc (i.e. trackres)
        %   column8 = x displacement based on gaussian fitting
        %   column9 = y displacement based on gaussian fitting
        %   column10= linear displacement based on gaussian fitting
        %   column11 = x displacement rmse of stdev
        %   column12 = y displacement rmse of stdev
        %   column13= linear displacement rmse of stdev w/ err propag
        
        tracksdisp(i,1)=xdisp;
        tracksdisp(i,2)=ydisp;
        tracksdisp(i,3)=lindisp;
        tracksdisp(i,4)=(mean(xdispmean{i}(:,2))-mean(xdispmean{i}(:,1)));
        tracksdisp(i,5)=(mean(ydispmean{i}(:,2))-mean(ydispmean{i}(:,1)));
        tracksdisp(i,6)=sqrt((mean(xdispmean{i}(:,2))-mean(xdispmean{i}(:,1)))^2+(mean(ydispmean{i}(:,2))-mean(ydispmean{i}(:,1)))^2);
        tracksdisp(i,7)=trackcalc{i}(1,8);
        tracksdisp(i,8)=dx;
        tracksdisp(i,9)=dy;
        tracksdisp(i,10)=dl;
        tracksdisp(i,11)=sx;
        tracksdisp(i,12)=sy;
        tracksdisp(i,13)=sl;
    end
end

%% write the errors and means into a matrix from the cell array
    % matrix for new cell array values
    % col 1,6  init,ext x position mean
    % col 2,7  init,ext y position mean
    % col 3,8  init,ext xx variance
    % col 4,9  init,ext yy variance
    % col 5,10 init,ext xy RMSE calculated
    
    % define a parameter for minimum motion for RMSE
    stuck=0.75;
    
    % enter the (error) matrix
    trackserror=[length(trackcalc) 14];
    
    %   set the conversion factor to go from pixels to microns
        % conv=-0.645;  %10X
        % conv=-1.625; %4x
         conv=-conv;
    % loop for extracting cell array values
    
    clear n;
    for n=1:length(trackcalc)
        % determine if cell is empty, continue if true
        if isempty(trackcalc{n,1});
            continue
        else
        % initial position
        trackserror(n,1)=trackfit{n,3}(1,1);    %meanx
        trackserror(n,2)=trackfit{n,3}(1,2);    %meany
        trackserror(n,3)=trackfit{n,4}(1,1); %variance XX
        trackserror(n,4)=trackfit{n,4}(2,2); %variance YY
        trackserror(n,5)=sqrt(...
            trackfit{n,4}(1,1) +...
            trackfit{n,4}(2,2));             %RMSE XY
        
        %extended position
        trackserror(n,6)=trackfit{n,6}(1,1);    %meanx
        trackserror(n,7)=trackfit{n,6}(1,2);    %meany
        trackserror(n,8)=trackfit{n,7}(1,1); %variance XX
        trackserror(n,9)=trackfit{n,7}(2,2); %variance YY    
        trackserror(n,10)=sqrt(...           %RMSE XY
            trackfit{n,7}(1,1) +...
            trackfit{n,7}(2,2));             
        
        % error propagation data
        trackserror(n,11)=sqrt(...           % X displacement error
            trackserror(n,3)^2 +...
            trackserror(n,8)^2);            
        trackserror(n,12)=sqrt(...           % Y displacement error
            trackserror(n,4)^2 +...
            trackserror(n,9)^2);            
        trackserror(n,13)=sqrt(...           % Lin displacement error
            trackserror(n,5)^2 +...
            trackserror(n,10)^2);;            
        
        %particle ID
        trackserror(n,14)=trackcalc{n,1}(1,1); %particle ID
        end
    end
    
    %   exclude obvious outliers
    show=[];
    for o=1:13  
        show{o, 1}=abs(trackserror(:,o))<2;
    end
    
    figure
    hold on
    histogram(trackserror(show{5},5).*conv^2,...
        'Normalization','probability','FaceColor',[0.7 0.7 0.7],...
        'BinWidth',.02)
    histogram(trackserror(show{10},10).*conv^2,...
        'Normalization','probability','FaceColor',[1 0 1],...
        'BinWidth',.02)
    histogram(trackserror(show{13},13).*conv^2,...
        'Normalization','probability','FaceColor',[0 0 1],...
        'BinWidth',.02)
    legend('initial position','extended position','displacement')
    ylabel('probability')
    xlabel('variance,microns^2')
    title(['RMSE Histogram fov ' num2str(fovn) ' ID ' id])
    hold off
    
     savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackAnalysis_RMSE-ini-ext-disp-Hist_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackAnalysis_RMSE-ini-ext-disp-Hist_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);
            
%% Determine the uniformity of the intial position

% define the parameter for symmetry
    symmetric=1.16;

for i=1:size(trackfit)
    if size(trackfit{i,4}) <= 1
                continue
    else
        %calculate the eigenvalue for the covariance matrix for initial
        %position infomation
        mx=max(eig(trackfit{i,4}));
        mn=min(eig(trackfit{i,4}));
        
        % save the key parameters
        % column 1 is the maximum eigenvalue of covariance matrix
        % column 2 is the minimum eigenvalue of the covariance matrix
        % column 3 is the uniformity parameter root of the ratio of  max to  min eigenvalue
        %           if column 3 (uniformity) is 1 it is perfectly uniform
        %           (i.e. max = min and max/min = 1)
        
        uniformity(i,1)=mx;     % max eigenvalue
        uniformity(i,2)=mn;     % min eigenvalue
        uniformity(i,3)=(mx/mn)^0.5; % this value <= the symmetry value is acceptable
    end
end

figure
    hold on
    % plot a histogram of the uniformity parameter
    histogram(uniformity(:,3),...%.*conv^2,...
        'Normalization','probability','FaceColor',[0.7 0.7 0.7],...
        'BinWidth',.05)
    legend('initial position')
    ylabel('probability')
    xlabel('symmetry parameter, microns^2')
    title(['Uniformity Parameter fov ' num2str(fovn) ' ID ' id])
    hold off
    
     savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackAnalysis_Uniformity-Hist_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackAnalysis_Uniformity-Hist_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);
%% Validated Tracks - save tracks with sufficient displacement 
% choose to comment/uncomment whether to analyze x or linear displacement

% %(i.e. legitimate tracks) to a new array, using x-displacement, single
% %threshold
%   toolegit=8; 
%     tracksreslegit={};
%         legit=find(abs(tracksdisp(:,4))>toolegit); 
%         tracksreslegit=tracksres(legit);
%      tracksresrejects={};
%           rejected=find(abs(tracksdisp(:,4))<=toolegit); 
%         tracksresrejects=tracksres(rejected);
        
 %(i.e. legitimate tracks) to a new array, using linear-displacement and
 %two thresholds, upper threshold is toolegit and lower is toquit
 
%  %toolegit is 7.5 to toquit 2 for 1pN
%     freakish=13.5;
%     toolegit=7.5;
%     toquit=2;
 

%  %toolegit 3 to toquit 1 for 3.6pN
%     toolegit=4;
%     toquit=1.5;

%  %toolegit 2 to toquit 0.5 for 5.4pN
%     toolegit=3;
%     toquit=1;

    tracksreslegit={};
        legitL=find(abs(tracksdisp(:,10))>=toolegit); 
        legitU=find(abs(tracksdisp(:,10))<freakish);
        legit=intersect(legitL,legitU);
        tracksreslegit=tracksreslonggated(legit);
    
    tracksresmediocre={};

        mehU=find(abs(tracksdisp(:,10)) < toolegit);      
        mehL=find(abs(tracksdisp(:,10)) > toquit);
        meh=intersect(mehU,mehL);
        tracksresmediocre=tracksreslonggated(meh);
               
    tracksresrejects={};
        rejected=find(abs(tracksdisp(:,10))<=toquit); 
        tracksresrejects=tracksreslonggated(rejected);      
        
    tracksresuniform={};
        uniform=find(uniformity(:,3)<=symmetric);
        tracksresuniform=tracksreslonggated(uniform);
        irregular=find(uniformity(:,3)>symmetric);
        
    tracksresmobile={};
        speedy=find(trackserror(:,10)>=stuck);
        quagmire=find(abs(trackserror(:,10))>stuck);
        tracksresmobile=tracksreslonggated(speedy);

 %% Exclusion of particles based on start time
%     
%     %define cell array for particles starting at 0 (the gate)
%     tracksresgated={};
%     tracksreslegitgated={};  
%     tracksresrejectsgated={};
%     tracksresmediocregated={};
%     
%     %define cell array for particles that arrive after the gate
%     tracksreslate={};
%     tracksreslegitlate={};   
%     tracksresrejectslate={};
%     tracksresmediocrelate={};
%     
%     % iterate through the tracks and write particles starting at the gate
%     % into tracksresgated cell array for all tracks
%     for i = 1: length(tracksres)
%        if tracksres{i}(1,7)==0
%            tracksresgated{i,1}=tracksres{i};    %write gated tracks 
%        else
%            tracksreslate{i,1}=tracksres{i};     %write late starting tracks
%        end
%     end
%     
%     % iterate through the tracks and write particles starting at the gate
%     % into tracksreslonggated cell array for all tracks
%     for i = 1: length(tracksreslong)
%        if tracksreslong{i}(1,7)==0
%            tracksreslonggated{i,1}=tracksreslong{i};    %write gated tracks 
%        else
%            tracksreslonglate{i,1}=tracksreslong{i};     %write late starting tracks
%        end
%     end
%     
%     % iterate through the tracks and write particles starting at the gate
%     % into tracksresgated cell array for legitimate tracks
%     for i = 1: length(tracksreslegit)
%        if tracksreslegit{i}(1,7)==0
%            tracksreslegitgated{i,1}=tracksreslegit{i};    %write legitimate gated tracks 
%        else
%            tracksreslegitlate{i,1}=tracksreslegit{i};     %write late legitimate starting tracks
%        end
%     end
%     
%     % iterate through the tracks and write particles starting at the gate
%     % into tracksresgated cell array for rejected tracks
%     for i = 1: length(tracksresrejects)
%        if tracksresrejects{i}(1,7)==0
%            tracksresrejectsgated{i,1}=tracksresrejects{i};    %write rejected gated tracks 
%        else
%            tracksresrejectslate{i,1}=tracksresrejects{i};     %write late rejected starting tracks
%        end
%     end
%     
%     % iterate through the tracks and write particles starting at the gate
%     % into tracksresgated cell array for mediocre tracks
%     for i = 1: length(tracksresmediocre)
%        if tracksresmediocre{i}(1,7)==0
%            tracksresmediocregated{i,1}=tracksresmediocre{i};    %write mediocre gated tracks 
%        else
%            tracksresmediocrelate{i,1}=tracksresmediocre{i};     %write mediocre legitimate starting tracks
%        end
%     end
%     

%% Save the variables that were calculated in this section
% legit tracks, gated tracks, long tracks, displacement tracks
save(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
    ['TrackAnalysis_fov_' num2str(fovn) '_v' num2str(ver)]),...
    'tracksdisp','trackserror','trackcalc','trackfit','dex','din',...
    'toolegit','toquit','tracksres', 'tracksresgated','tracksreslate',...
    'tracksreslong','tracksreslonglate','tracksreslonggated', ...
    'tracksreslegit', 'legit',...
    'tracksresrejects',...
    'rejected','tracksresmediocre',...
    'meh','uniformity','stuck',...
    'tracksresuniform','uniform','irregular','tracksresmobile',...
    'speedy','quagmire');

%'tracksreslegitgated','tracksreslegitlate',...
%'tracksresrejectslate','tracksresrejectsgated',...
%'tracksresmediocregated','tracksresmediocrelate'
%% Plot new analysis of displacement in pixels

%   set the number of bins for the histogram
        l=round(sqrt(length(tracksdisp)));  
        biwi=0.3225; %bin width

%   exclude obvious outliers
        showx=abs(tracksdisp(:,8))<20;
        showy=abs(tracksdisp(:,9))<20;
        showl=abs(tracksdisp(:,10))<20;
%   specify which type of averaging is used
    class='gausian';

    % Plot histogram overlay in pixels of all displacements
    figure('Name',['X Y Linear Displacement histogram / ' class]);
    hold on
    hx=histogram(tracksdisp(showx,8).*-1,'BinWidth', abs(biwi/conv), 'Normalization','count','FaceColor','r');
    hy=histogram(tracksdisp(showy,9).*-1,'BinWidth', abs(biwi/conv), 'Normalization','count','FaceColor','b');
    hl=histogram(tracksdisp(showl,10).*1,'BinWidth', abs(biwi/conv), 'Normalization','count','FaceColor',[0 0 0]);
        
            title(['Track Analysis X, Y, Linear Displacement Histogram Stack ' name ' fov' num2str(fovn)]);
            legend('X' , 'Y', 'Linear');            
            xlabel('displacement,  pixels');
            ylabel('events');
             xlim([-15 15])
    hold off
    
    % (un)comment these lines to save the x displacement file individually
        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_XYLinpixel-DisplacementHistogram_' name '_fov' num2str(fovn) '.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_XYLinpixel-DisplacementHistogram_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);

    % Plot histogram of the x displacement of all tracks
    figure('Name',['X Displacement histogram / ' class]);

        hx=histogram(tracksdisp(showx,8).*conv,'BinWidth', biwi, 'Normalization','count','FaceColor','r');
            title(['Track Analysis X Displacement Histogram Stack ' name ' fov' num2str(fovn)]);
            legend(['x-displacement']);            
            xlabel('displacement,  \mum');
            ylabel('events');
             xlim([-2 12])


    % (un)comment these lines to save the x displacement file individually
        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_X-DisplacementHistogram_' name '_fov' num2str(fovn) '.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_X-DisplacementHistogram_' name '_fov' num2str(fovn) '.png']),'-dpng');
        
    close(gcf);

    % % %(un)comment this line to plot x and y disp individuallt
    %     hold on;    %hold to plot y displacements on the same bar graph

% Plot histogram of the maximum x displacement of all tracks
    figure('Name',['X Displacement histogram / ' class]);

        hy=histogram(tracksdisp(showy,9).*conv, 'BinWidth', biwi,'Normalization','count','FaceColor','b');
    
            title(['Track Analysis Y Displacement Histogram Stack ' name ' fov' num2str(fovn)]);
            legend(['y-displacement']);           
            xlabel('displacement, \mum');
            ylabel('events');
             xlim([-2 12])

    % (un)comment these lines to save the x displacement file individually
        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_Y-DisplacementHistogram_' name '_fov' num2str(fovn) '.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_Y-DisplacementHistogram_' name '_fov' num2str(fovn) '.png']),'-dpng');

    %         savefig(fullfile(basepath, 'Bead_tracking',['TrackAnalysis_X&Y-DisplacementHistOverlay_fov' num2str(fovn) '.fig']));
    %         print('-f',fullfile(basepath, 'Bead_tracking',['TrackAnalysis_X&Y-DisplacementOverlay_fov' num2str(fovn) '.png']),'-dpng');

    close(gcf);
    
% Plot histogram of the linear displacement of all tracks
    figure('Name',['Linear Displacement histogram / ' class]);

        hl=histogram(tracksdisp(showl,10).*-conv, 'BinWidth', biwi,'Normalization','count','FaceColor',[0 0 0]);

            title(['Track Analysis Linear Displacement Histogram Stack ' name ' fov' num2str(fovn)])
            legend(['linear displacement']);           
            xlabel('displacement, \mum');
            ylabel('events');
             xlim([-2 12])

    % (un)comment these lines to save the x displacement file individually
        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_Lin-DisplacementHistogram_' name '_fov' num2str(fovn) '.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_Lin-DisplacementHistogram_' name '_fov' num2str(fovn) '.png']),'-dpng');

    %         savefig(fullfile(basepath, 'Bead_tracking',['TrackAnalysis_X&Y-DisplacementHistOverlay_fov' num2str(fovn) '.fig']));
    %         print('-f',fullfile(basepath, 'Bead_tracking',['TrackAnalysis_X&Y-DisplacementOverlay_fov' num2str(fovn) '.png']),'-dpng');

    close(gcf);

% Plot as ovelay x and y displacements
figure('Name',[' Displacement histogram stack / ' class]);

    %plot x displacements
    subplot(3,1,1)
        hx=histogram(tracksdisp(showx,8).*conv,'BinWidth', biwi,'Normalization','probability','FaceColor','r');
             xlim([-2 12])
                %title(['Track Analysis X&Y Displacement Histogram Stack ' name ' fov' num2str(fovn)])
                legend(['\itx\rm-displacement'])            
                xlabel('\itx\rm-displacement, \mum');
                ylabel('frequency'); 
   

   
    % plot y displacements
    subplot(3,1,2)
        hy=histogram(tracksdisp(showy,9).*conv, 'BinWidth', biwi,'Normalization','probability','FaceColor','b');
             xlim([-2 12])
                legend(['\ity\rm-displacement'])            
                xlabel('\ity\rm-displacement, \mum');
                ylabel('frequency'); 
   
                
    % plot linaer displacements
    subplot(3,1,3)
        hl=histogram(tracksdisp(showl,10).*-conv, 'BinWidth', biwi,'Normalization','probability','FaceColor',[0 0 0]);
             xlim([-2 12])
                legend(['linear-displacement'])            
                xlabel('displacement, \mum');
                ylabel('frequency'); 
 

        savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
            ['TrackAnalysis_X&Y-DisplaceHistStack_' name '_fov' num2str(fovn) '.fig']));
        print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
        ['TrackAnalysis_X&Y-DisplacementStack_' name '_fov' num2str(fovn) '.png']),'-dpng');
     
    close(gcf);

%% plot the individual distributions
figure
hold on
histogram(tracksdisp(legit,10).*-conv,...
    'BinWidth', abs(biwi*conv),'Normalization','count','FaceColor','r');
histogram(tracksdisp(meh,10).*-conv,...
    'BinWidth', abs(biwi*conv),'Normalization','count','FaceColor',[0 0.4 0.4 ]);
histogram(tracksdisp(rejected,10).*-conv,...
    'BinWidth', abs(biwi*conv),'Normalization','count','FaceColor',[0 0.4 0]);
xlabel('linear displacement')
ylabel('count')
legend('legitimate','mediocre','rejected')
xlim([0 10])    
hold off
     savefig(fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackAnalysis_LMR-DisplaceHistIndiv_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackAnalysis_LMR-DisplacementIndiv_' name '_fov' num2str(fovn) '.png']),'-dpng');
     
    close(gcf)
    
    %% fit legit particles to a normal distribytion
% 
%     figure('Name',['Fit PDF to Linear Displacement histogram / ' class]);
%     hold on
%     clear Data
%     k=1;
%     sign=1;
%     direct=' ';
%     labels={};
%     
%     for i=8:10
%         if i==4,8;
%             direct='X';
%             sign=1;
%         elseif i==5,9;
%             direct='Y';
%             sign=1;
%         elseif i==6,10;
%             direct='Linear';
%             sign=-1;
%         elseif i==4,5,6;
%             class='1D';
%         elseif i==8,9,10;
%             class='gaussian';
%         else i==7
%             continue
%         end
%             
%         for j=1:3
%             if j==1;
%                 set=legit;
%                 setname='Legitimate';
%             elseif j==2;
%                 set=meh;
%                 setname='Meidocre';
%             else j==3;
%                 set=rejected;
%                 setname='Rejects';
%             end
%                 
%             Data=tracksdisp(set,i).*conv.*sign;
%         
%             h=histogram(Data, 'NumBins', l,'Normalization','pdf');
% 
%                 title(['Track Analysis ' direct '-Displacement Histogram Stack ' setname ' fov' num2str(fovn)]);
%                 labels{i}=[direct ' ' setname ' ' class ];
%                 xlabel('displacement, \mum');
%                 ylabel('density');
%                 xlim([-2 10]);
%            
%             if j==2;
%                 continue
%             elseif j==3;
%                 continue
%             else j==1;
%             % Create grid where function will be computed
%             XLim = [-2 10];
%             XLim = XLim + [-1 1] * 0.01 * diff(XLim);
%             XGrid = linspace(XLim(1),XLim(2),100);
% 
%             pd1 = fitdist(Data, 'normal');
%             YPlot = pdf(pd1,XGrid);
%             plot(XGrid,YPlot,'LineWidth',2);
% 
%             clear fitLMR
%             fitLMR={};
% 
%             fitLMR{k,1}=[setname '_' direct '_' class '_microns'];  %name for the data
%             fitLMR{k,2}=Data;             % save the dat that was fit
%             fitLMR{k,3}=YPlot;            % save the y data for the fit
%             fitLMR{k,4}=XGrid;            % save the x data for the fit
%             fitLMR{k,5}=pd1;              % save the prob density function fit
%             fitLMR{k,6}(1,1)=mean(pd1);   % obtain the coefficients mu and sigma
%             fitLMR{k,6}(1,2)=std(pd1);
%             fitLMR{k,7}=paramci(pd1);     % obtain the confidence itnervals
%             k=k+1;
%             end
%         end
%     end
%     legend({'Legitimate';['Fit mu=' num2str(fitLMR{1,6}(1,1)) ' sig=' num2str(fitLMR{1,6}(1,2)) ]; ...
%         'Mediocre';'Rejected'})
%     hold off
%     clear k
%         
%         %save the figure with color coded bins, fit to accepted tracks
%         savefig(fullfile(basepath, 'Bead_tracking',['TrackAnalysis_SegregatedHistogram_' name '_fov' num2str(fovn) '.fig']));
%         print('-f',fullfile(basepath, 'Bead_tracking',['TrackAnalysis_SegregatedHistogram_' name '_fov' num2str(fovn) '.png']),'-dpng');
%      
%     close(gcf);
%     
%     %save the FitLMR variable
%     save(fullfile(basepath,'Bead_tracking',['TrackAnalysis_Fits_fov' num2str(fovn)])...
%         ,'fitLMR');

%% end
%version1 - simply segregated particles above and below a threshold
%version2 - 2/2018 use two thresholds to bin the particles further
%version3 - 2/28/2018 plot histograms in microns
%version4 - 3/2/2018 plot normalized histograms
%version5 - 5/11/2018 made the thresholds a variable for input
%version6 - 5/31/2018 explicitly calculate the error based on gaussian,
%                       use version to create a TrackAnalysis_Sorting dir
%version7 - 7/5/2018 include a calculatin of uniformity and mobility
%                    based on the initial position data, sort based on this