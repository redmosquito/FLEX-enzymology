%% count particles from tracks trajectory file containing frame numbers
% input is a string corresponding to the particles to be tracked
% this code will be ran after tracking all the particles. Code needs to be
% run from the root directory that initially contained the image file
%
% Drew Drabek
% Blacklow and Loparo Labs
% Harvard Medical School, BCMP Dept
%
% 2018/02/13 AD - Updated with capacity to create three bins when sorting
% tracks for counting

function particlecounter_complete(fovn,ver,fr)
%% PC Initialize the path and indicate the version of tracks to be analyzed
%  %Document the working directory for the images
%     basepath = [cd,'\'];
%     
% % initialize the script with a few constants
% 
% % % Enter an index for saving additional tracked movies (for repeated
% % % analysis of the same set of trajectories)
% %     ver=1;
% % 
% % % Enter a frame rate for plotting in terms of time signature
% %     fr=1; %frames per second
%     
% % Define the conversion from pixels to microns
%     conv=-0.645;
%         
% %Create a folder for the analysis
%     mkdir(basepath,['Track_analysis_v' num2str(ver)]);
%     
%     load(fullfile(basepath, 'Bead_tracking\', ...
%         ['TrackAnalysis_fov' num2str(fovn) '.mat']));
%     
%     load(fullfile(basepath, 'Bead_tracking\', ...
%         ['TrackAnalysis_params_fov' num2str(fovn) '.mat']));
%     
%     load(fullfile(basepath, 'Bead_tracking\res_files\', ...
%         ['fov' num2str(fovn) '_pretrackparams.mat']));
%     
%% MAC Initialize the path and indicate the version of tracks to be analyzed
 %Document the working directory for the images
    basepath = [cd ''];
    
% initialize the script with a few constants

% % Enter an index for saving additional tracked movies (for repeated
% % analysis of the same set of trajectories)
%     ver=1;
% 
% % Enter a frame rate for plotting in terms of time signature
%     fr=1; %frames per second
    
% Define the conversion from pixels to microns
    conv=-0.645;
        
%Create a folder for the analysis
    mkdir(basepath,['Track_analysis_v' num2str(ver)]);
    
    load(fullfile(basepath, ['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)], ...
        ['TrackAnalysis_fov_' num2str(fovn) '_v' num2str(ver) '.mat']));
        basepath = [cd ''];

    load(fullfile(basepath, ['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)], ...
        ['TrackAnalysis_params_fov_' num2str(fovn) '_v' num2str(ver) '.mat']));
        basepath = [cd ''];

    load(fullfile(basepath, 'Bead_tracking/res_files/', ...
        ['fov' num2str(fovn) '_pretrackparams.mat']));
        basepath = [cd ''];

%% Initialize the script by choosing an input track set and name (manual)
% enter the name of the input tracks cell array
% enter a name for a handle on files and figures for this counting
    
    % plotting of particles based on mean square deviation
    %    input = tracksresmoving;
    %     name = 'MovingTracks';
    %    input = tracksresstill;
    %     name = 'StillTracks';

%     %unsorted trajectories
%     input = tracksres;
%     name = 'AllTracks_unsorted';
%     j=0; k=uns;

%      %tracks  to calculate displacements, greater than min time (long)
%     input = tracksreslonggated;
%     name = 'AllTracks';
%     j=1; k='all';

%     %tracks that are long, start at time zero, move a key distance
%          input = tracksreslegitgated;
%          name = 'LegitTracks_7pt5Linear';
%          j=2; k='leg';

%      input = tracksresmehgated;
%      name = 'MediocreTracks_7pt5-2pLinear';
%      j=3; k='med';
     
%      input = tracksresrejectsgated;
%      name = 'RejectedTracks_2pt0Linear';
%      j=4; k='rej';
%      

%% Load the track analysis from the Track_sorting folder

%     trackanalysis=load(fullfile(basepath, ['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)], ...
%         ['TrackAnalysis_fov' num2str(fovn) '.mat']));
    trackanalysis=load(fullfile(basepath, ['Track_sorting_fov' num2str(fovn) '_v' num2str(ver)], ...
        ['TrackAnalysis_fov_' num2str(fovn) '_v' num2str(ver) '.mat']));
    
    tracknames=fieldnames(trackanalysis);
    tocountnames=strfind(tracknames,'tracksres','ForceCellOutput',true);
    nocount=strfind(tracknames,'late','ForceCellOutput',false);
    
    clear trackstocount
    trackstocount=cell(1,1);
    j=1;
    for i = 1:length(tracknames)
        if isempty(nocount{i})==0
            continue
        elseif  isempty(tocountnames{i}) == 0
            
            % assign the field ith name index into jth row 1st col
            trackstocount{j,1}=getfield(trackanalysis,tracknames{i,1});
            
            % assign the name of the ith field into jth row 2nd col
            trackstocount{j,2}=tracknames{i,1};
            
            % increase the index by 1 for the assignment
            j=j+1;
        else
                continue
                %[ 'trackstocount ' num2str(i) ' is empty'] %command line output
        end
   end 
    clear j
    
%% assign known variables to handles for plotting

clear i
clear k
clear remaining

%define the cell array to summarize the data
remainingsummary=cell(length(trackstocount),4);

for k=1:length(trackstocount)
    input=trackstocount{k,1};
    name=trackstocount{k,2};
    
%% Perform counting
% total number of particles counted
    total=length(input);
    runningtotal=total;
% define a variable for the beads remaining
    remaining=[];
    remaining=[runningtotal 0];

% count the total frames for the trajectory
    lastframe=cellfun('length',input);
    numframes=max(lastframe(:,1));
    %numframes=60;  %uncomment to test with only 60 frames
    
%% use a for loop to subtract particle for each track ending
    for i=1:numframes
        numlost = find(lastframe(:,1) == i);
        runningtotal = runningtotal - length(numlost);
        remaining(i,1)=runningtotal;
        remaining(i,2)=i;
    end
    
    % calculate the fraction beads remaining
    remaining(:,3)=remaining(:,1)./max(remaining(:,1));
    
    % save the remaining data into a cell array
    remainingsummary{k,1}=remaining;
    remainingsummary{k,2}=name;
    remainingsummary{k,3}=input;
    %remainingsummary{k,2}=inputs{2,1}(1,k);
    if size(remainingsummary{k,3}) <= 1
            continue
    end
        
%% save parameters used for counting
    save(fullfile(basepath,['Track_analysis_v' num2str(ver)],...
        [name '_remaining_fov' num2str(fovn) '_v' num2str(ver) ]),...
        'remaining');

    csvwrite(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
        [name '_remaining_fov' num2str(fovn) '_v' num2str(ver) '.csv']),...
        remaining);
% clear remaining
%% Plot Results, Save figures

    % plot results and save a figure
    figure('Name',name)
        plot(remaining(:,2)./fr,remaining(:,1),'ko')
            ymin=0;
            ymax=max(remaining(:,1));
            ylim([ymin ymax])
        grid on
        xlabel('time, seconds');
        ylabel('number of particles');
   
    savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
        [name '_countVtime_fov' num2str(fovn) '_v' num2str(ver) '.fig']));
    print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
        [name '_countVtime_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');

    close(gcf);
 
    % Fraction Beads vs Time plot results and save a figure
        figure('Name',name)
            hold on
                     
            plot(remaining(:,2)./fr,remaining(:,3),'ko')
            
            grid on
            xlabel('time, seconds');
            ylabel('fraction of particles remaining');
                ymin=0;
                ymax=1;
                ylim([ymin ymax]);
        savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
            [name '_fractionVtime_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
        print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
            [name '_fractionVtime_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');

        close(gcf);
        
end

% save the summary of the bead loss trajectories (remaining)
    save(fullfile(basepath,['Track_analysis_v' num2str(ver)],...
        ['Track_analysis_remainingsummary_fov' num2str(fovn) '_v' num2str(ver) ]),...
        'remainingsummary');

%     csvwrite(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%         ['Track_analysis_remainingsummary_fov' num2str(fovn) '_v' num2str(ver) '.csv']),...
%         remainingsummary);
    
%% Definte key sets of trajectories to be plotted
    clear inputs
    inputs=cell(2,1);
    inputs{1,1}={...
        'tracksres',...
        'tracksreslonggated',...
        'tracksreslegit',...
        'tracksresmediocre',...
        'tracksresrejects'};
    inputs{2,1}={...
        'AllTracks_unsorted', ...
        'AllTracks', ...
        ['LegitTracks_' num2str(round(toolegit*-conv,1)) 'um'], ...
        ['MediocreTracks_' num2str(round(toolegit*-conv,1)) '-' num2str(round(toquit*-conv,1)) 'um'], ...
        ['RejectedTracks_' num2str(round(toquit*-conv,1)) 'um']};
%% overlay the key trajectories without fitting
% Fraction Beads vs Time plot results and save a figure
                     
    figure('Name','Overlay All, Legit, Mediocre, Rejected')
        hold on
            
        % define a color scheme for the plots   
            colors= {   'k.',...
                        'k-',...    
                        'r-',...    
                        'c-',...    
                        'b-',...    
                            };
        clear i
        clear j
        clear id
        delay=extend+extlen;
        % save the names from remainingsummary to their own cell array
        clear names
            names=cell(1);%define a cell array for saving names
            for i=1:length(remainingsummary)
                names{i,1}=remainingsummary{i,2};
            end
            
        % define the number of terminal frames to exclude
        j=5;
        
       for i=1:5
           if i==1
               continue
           
                         
           else
            clear id idname
            c=colors{i};
            id=inputs{1,1}(1,i);
            idname=id{1};
            %idname=idname.name
                
            toplot=find(strcmp(names,idname));
            
%            if  size(remainingsummary{toplot,3}) <= 1
%                 continue
%            end
           
            x=remainingsummary{toplot,1}(1:end-j,2)./fr;
            y=remainingsummary{toplot,1}(1:end-j,3);
           
            %plot(x(1:end-delay+1),y(delay:end),c)
            plot(x,y,c)
           end
        end
                xlabel('time, seconds');
                ylabel('fraction of particles remaining');
                ymin=0;
                ymax=1;
                ylim([ymin ymax]);
        legend('all','legitimate','mediocre','rejected')
        hold off

        savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
            ['Overlay_fractionVtime_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
        print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
            ['Overlay_fractionVtime_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');

        close(gcf);
%% overlay the key trajectories with 2phase exponential fitting
% % Fraction Beads vs Time plot results and save a figure
%                  startfit=extend+extlen;
%                  endfit=startfit+2000;
%     fits=cell(5,3);
%      % define a color scheme for the plots   
%             colors= {   'k',...
%                         'k',...    
%                         'r',...    
%                         'c',...    
%                         'b',...    
%                             };
%     figure('Name','Overlay All, Legit, Mediocre, Rejected with Fit')
%         hold on
%             
%         clear i
%         clear j
%         clear id
%         delay=extend+extlen;
%         % save the names from remainingsummary to their own cell array
%         clear names
%             names=cell(1);%define a cell array for saving names
%             for i=1:length(remainingsummary)
%                 names{i,1}=remainingsummary{i,2};
%             end
%             
%         % define the number of terminal frames to exclude
%         j=5;
%         
%         for i=2:5
%            
%             clear id idname
%             c=colors{i};
%             id=inputs{1,1}(1,i);
%             idname=id{1};
%             %idname=idname.name
%                 
%             toplot=find(strcmp(names,idname));
%             x=remainingsummary{toplot,1}(1:end-j,2)./fr;
%             y=remainingsummary{toplot,1}(1:end-j,3);
%             fx=remainingsummary{toplot,1}(startfit:endfit,2)./fr;
%             fy=remainingsummary{toplot,1}(extend+extlen:end-j,3);
%             
% %             f=fit(fx,fy,'exp2')         % fit 2 phase exponential decay
%             f=fit(fx,fy,'exp2','Exclude',fx < extend+extlen );         % fit 2 phase exponential decay
%             fits{i,1}=f;                % save the fit to a cell
%             
%             fits{i,2}=coeffvalues(f);    % save the coefficients to a cell
%                                         % a N01 b lambda1 / c N02 d lambda2
%             fits{i,3}=confint(f);        % save the confidence intervals
%              
%             if i==1
%                 continue
%             else
%             p=plot(f,fx,fy);
%             p(1).Color=c;
% %             p(1).Line='--';
%             p(2).Color=c;
%             %p(2).marker='d'
% 
%             end      
%         end
%                 xlabel('time, seconds');
%                 ylabel('fraction of particles remaining');
%                 ymin=0;
%                 ymax=1;
%                 ylim([ymin ymax]);
%         legend('all','all fit','legit','legit fit','mediocre','med fit','rejected','rej fit')
%         hold off
% 
%         savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%             ['Overlay_fractionVtime__Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
%         print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%             ['Overlay_fractionVtime__Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');
% 
%         close(gcf);
%         
%         % save the fits
%          save(fullfile(basepath,['Track_analysis_v' num2str(ver)],...
%         ['TrackAnalysis_Fits_fov' num2str(fovn) '.mat']),...
%         'fits');
%     
    %% Plot a pie chart of the total particles in each sorted bin
    labels={'FULL',...
        'PART',...
        'IMMOB'};
    slices=[length(remainingsummary{5,3})...
        length(remainingsummary{6,3}) ...
        length(remainingsummary{7,3})];
     figure
    % hold on
    hpie=pie(slices,[1 0 0]);
    hpie(1).FaceColor='red';
    hpie(2).FontSize=14;
    hpie(3).FaceColor='blue';
    hpie(4).FontSize=14;
    hpie(5).FaceColor='cyan';
    hpie(6).FontSize=14;
    legend(labels,'location','best')
    %hold off
    
        savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
            ['Piechart_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
        print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
            ['Piechart_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');

        close(gcf);
    
%% Fit Data to a single-molecule kinetic model

% 
% % d[P]/dt=(([E]*kcat)/(KD+[E])*(C0-P))
% E=10E-6;                         %[E] in M
% k2=0.02;                         % kcat in s^-1
% KD=10E-6;                        %dissociation constant in M
% C=100; %re-defined in the loop below %number of monovalently tethered particles
% 
% %E  = [E] in M - constant approx 0.02s-1
% %KD = KD in M - constant approx 5-15uM
% %k2 = kcat in s^-1
% %C  = C0, number of particles in first frame
% %P  = number of product molecules at time t 
% 
% 
% % perform nonlinear fit to a single molecule kinetics equation
%  fits=cell(5,3);
%      % define a color scheme for the plots   
%             colors= {   'k',...
%                         'k',...    
%                         'r',...    
%                         'c',...    
%                         'b',...    
%                             };
%     figure('Name','Overlay All, Legit, Mediocre, Rejected with Nonlinear Fit')
%         hold on
%             
%         clear i
%         clear j
%         clear id
%         delay=extend+extlen;
%         % save the names from remainingsummary to their own cell array
%         clear names
%             names=cell(1);%define a cell array for saving names
%             for i=1:length(remainingsummary)
%                 names{i,1}=remainingsummary{i,2};
%             end
%             
%         % define the number of terminal frames to exclude
%         j=5;
% i=2;
%     clear id idname
%     c=colors{i};
%     id=inputs{1,1}(1,i);
%     idname=id{1};
%     
%     toplot=find(strcmp(names,idname));
%     C=remainingsummary{toplot,1}(1,1)
%     x=remainingsummary{toplot,1}(1:end-j,2)./fr;
%     y=remainingsummary{toplot,1}(1:end-j,3);
%     fx=remainingsummary{toplot,1}(extend+extlen:end-j,2)./fr;
%     fy=remainingsummary{toplot,1}(extend+extlen:end-j,3);
%     
%     fsm =@(b,P) ((b(2)*b(1))/(b(3)+b(1)))*(b(4)-P);
%     clear beta
%     beta=[ E; k2; KD; C]
%     bfit = nlinfit(fx,fy,fsm,beta)
% 
% 
% % example from the web
% % f = @(b,x) 3.379 - (b(1)*x.^2./(b(2)+x))
% % bfit = nlinfit(xdata,ydata,f,b0)
% for i=2:5
%            
%             clear id idname
%             c=colors{i};
%             id=inputs{1,1}(1,i);
%             idname=id{1};
%             %idname=idname.name
%                 
%             toplot=find(strcmp(names,idname));
%             x=remainingsummary{toplot,1}(1:end-j,2)./fr;
%             y=remainingsummary{toplot,1}(1:end-j,3);
%             fx=remainingsummary{toplot,1}(extend+extlen:end-j,2)./fr;
%             fy=remainingsummary{toplot,1}(extend+extlen:end-j,3);
%             
% %             f=fit(fx,fy,'exp2')         % fit 2 phase exponential decay
%             f=fit(fx,fy,'exp2','Exclude',fx < extend+extlen )         % fit 2 phase exponential decay
%             fits{i,1}=f;                % save the fit to a cell
%             
%             fits{i,2}=coeffvalues(f)    % save the coefficients to a cell
%                                         % a N01 b lambda1 / c N02 d lambda2
%             fits{i,3}=confint(f)        % save the confidence intervals
%              
%             if i==1
%                 continue
%             else
%             p=plot(f,fx,fy)
%             p(1).Color=c;
% %             p(1).Line='--';
%             p(2).Color=c;
%             %p(2).marker='d'
% 
%             end      
%         end
%                 xlabel('time, seconds');
%                 ylabel('fraction of particles remaining');
%                 ymin=0;
%                 ymax=1;
%                 ylim([ymin ymax]);
%         legend('all','all fit','legit','legit fit','mediocre','med fit','rejected','rej fit')
%         hold off

%% load coefficients to be plotted into variables
%        General model Exp2:
%     f(x) = a*exp(b*x) + c*exp(d*x)
%   a and c are the initial quantity
%   b and d are decay rate

%     % load fast decay constants (lambda) into an array
%     fastlam=[fits{2,2}(1,2) fits{3,2}(1,2) fits{4,2}(1,2) fits{5,2}(1,2)].*-1; %decay rate, s-1
%     fasttau=1./fastlam; %time constant, s
%     fasthalf=0.6932./fastlam; %half life, s
%     fastnorm=abs(fits{2,2}(1,2));
%     
%     % load slow decay constants (lambda) into an array
%     slowlam=[fits{2,2}(1,4) fits{3,2}(1,4) fits{4,2}(1,4) fits{5,2}(1,4)].*-1; %decay rate, s-1
%     slowtau=1./slowlam; %time constant, s
%     slowhalf=0.6932./slowlam; %half life, s
%     slownorm=abs(fits{2,2}(1,4));
%     
%% plot decay rate, also called lambda and K
%        General model Exp2:
%     f(x) = a*exp(b*x) + c*exp(d*x)
% decay constant is b and d, units are s-1

%     label='decay constant, 1/s';
%     fast=fastlam;
%     slow=slowlam;
%        
%     % plot the rates for the fast and slow
%     figure; hold on
%         subplot(3,2,1)
%             bar(fast,'g')
%                 ylabel(label);
%                 cols={'all';'leg';'med';'rej'};
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 legend('fast');
%         subplot(3,2,3)
%             bar(slow,'y')
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylabel(label);
%                 legend('slow');
% 
%         subplot(3,2,5)
%             clear i
%             for i=1:length(fast)
%                 pl(i,1)=fast(i);
%                 pl(i,2)=slow(i);
%             end
%             
%             bar([pl(1,1) pl(1,2); pl(2,1) pl(2,2); pl(3,1) pl(3,2); pl(4,1) pl(4,2)],'stacked')
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylabel(label);
%                 legend('fast','slow')
%         
%         % plot the normalized rates for fast and slow
%          subplot(3,2,2)
%             bar([fastlam./fastnorm]-1,'g')
%                 ylabel(label);
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylim([-0.5 0.5])
%         subplot(3,2,4)
%             bar([slowlam./slownorm]-1,'y')
%                 ylabel(label);
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylim([-0.25 0.25])
%         hold off
%         
%         savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%             ['Overlay_DecayRates_K_Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
%         print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%             ['Overlay_DecayRates_K_Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');
% 
%         close(gcf);
% 
% %% plot exponential decay constant, also called mean lifetime, tau
% %        General model Exp2:
% %     f(x) = a*exp(b*x) + c*exp(d*x)
% %   mean lifetime is -1/b
% 
%     label='mean lifetime, s';
%     fast=fasttau;
%     slow=slowtau;
%     
%      % plot the rates for the fast and slow
%     figure; hold on
%         subplot(3,2,1)
%             bar(fast,'g')
%                 ylabel(label);
%                 cols={'all';'leg';'med';'rej'};
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 legend('fast');
%         subplot(3,2,3)
%             bar(slow,'y')
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylabel(label);
%                 legend('slow');
% 
%         subplot(3,2,5)
%             clear i
%             for i=1:length(fast)
%                 pl(i,1)=fast(i);
%                 pl(i,2)=slow(i);
%             end
%             
%             bar([pl(1,1) pl(1,2); pl(2,1) pl(2,2); pl(3,1) pl(3,2); pl(4,1) pl(4,2)],'stacked')
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylabel(label);
%                 legend('fast','slow')
%         
%         % plot the normalized rates for fast and slow
%          subplot(3,2,2)
%             bar([fasttau./(1./fastnorm)]-1,'g')
%                 ylabel(label);
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylim([-0.75 0.75])
%         subplot(3,2,4)
%             bar([slowtau./(1./slownorm)]-1,'y')
%                 ylabel(label);
%                 set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                 ylim([-0.5 0.5])
%         hold off
%         
% savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%             ['Overlay_DecayConstant_T_Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
%         print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%             ['Overlay_DecayConstant_T_Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');
% 
%         close(gcf);
%% plot half life     

%         label='half-life, s';
%         fast=fasthalf;
%         slow=slowhalf;
% 
%              % plot the rates for the fast and slow
%         figure; hold on
%             subplot(3,2,1)
%                 bar(fast,'g')
%                     ylabel(label);
%                     cols={'all';'leg';'med';'rej'};
%                     set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                     legend('fast');
%             subplot(3,2,3)
%                 bar(slow,'y')
%                     set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                     ylabel(label);
%                     legend('slow');
% 
%             subplot(3,2,5)
%                 clear i
%                 for i=1:length(fast)
%                     pl(i,1)=fast(i);
%                     pl(i,2)=slow(i);
%                 end
% 
%                 bar([pl(1,1) pl(1,2); pl(2,1) pl(2,2); pl(3,1) pl(3,2); pl(4,1) pl(4,2)],'stacked')
%                     set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                     ylabel(label);
%                     legend('fast','slow')
% 
%             % plot the normalized rates for fast and slow
%              subplot(3,2,2)
%                 bar([fasthalf./(0.6392./fastnorm)]-1,'g')
%                     ylabel(label);
%                     set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                     ylim([-0.75 0.75])
%             subplot(3,2,4)
%                 bar([slowhalf./(0.6392./slownorm)]-1,'y')
%                     ylabel(label);
%                     set(gca, 'XTickLabel',cols,'XTick',1:numel(cols))
%                     ylim([-0.5 0.5])
%             hold off
% 
%             savefig(fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%                 ['Overlay_HalfLife_Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.fig']));
%             print('-f',fullfile(basepath, ['Track_analysis_v' num2str(ver)],...
%                 ['Overlay_HalfLife_half_Fit2exp_fov' num2str(fovn) '_v' num2str(ver)  '.png']),'-dpng');
% 
%             close(gcf);
        
        %% pie chart
%     a=
%     m=
%     r=
%     labels={'>= 6.13uM','1.3 - 6.13um','<=1.3uM'};
%     pie([ all mediocre rejects],[1 0 0 ],labels)
%     
%% count crudely using number of objects per frame, unlinked (manual)
% this analysis is now included in the complete_tracking script 3/6/18
% 
% load(fullfile(basepath, 'Feature_finding/', ['MT_' num2str(fovn) '_Feat_Size_' num2str(featsize) '.mat']));
% 
% clear i indices
%     
%     % use a loop to save the number of accepted objects per frame
%     for i = 1:max(MT(:,6))
%         % assign the indices for the ith particle
%         indices=find(MT(:,6)==i);
%         
%         % for MT file save the number of objects per frame
%         frames{i}=[MT(indices(1):indices(end),7) MT(indices(1):indices(end),1) MT(indices(1):indices(end),2)];
%         
%         % for MT file save number of objects in each frame
%         % column 1 = time signature, ith frame
%         % column 2 = number of objects in the ith frame
%         framesct(i,1)=MT(indices(1),6);
%         framesct(i,2)= length(indices);
%         
%     end
% 
% % plot crude analysis (akin to 'ImageJ' analysis, but eliminated clusters)
%     figure;
%     fr=2;
%         hold on
%         T=['fov' num2str(fovn) ' Crude # Particles vs Time'];
%             xlabel('time, s');
%             ylabel('fraction of initial objects remaining');
%             %ymax=max(framesct(:,1)); %maximum for count
%             ylim([0 1]);
%             plot(framesct(:,2)./fr,framesct(:,1)./ymax)
%         hold off 
%             savefig(fullfile(basepath, 'Bead_tracking',['PlotFrac_raw_fov' num2str(fovn) '.fig']));
%             print('-f',fullfile(basepath, 'Bead_tracking',['PlotFracvs_raw_fov' num2str(fovn) '.png']),'-dpng');
%     close(gcf);
%     clear fr;
end

%% Appendix for comparative counting
% 
% % need to load the remaining vector for each set of sorted tracks into a
% % new variable
% rejectcount=remaining;
% count=remaining;
% legitcount=remaining;
% 
% % need to plot the variables on the same axis
% plot(count(:,2),count(:,3),rejectcount(:,2),rejectcount(:,3),legitcount(:,2),legitcount(:,3))
%     legend('AD08-049-1 all','AD08-049-1 rejects','AD08-049-1 legit')
%              grid on
%             xlabel('time, seconds');
%             ylabel('fraction of particles remaining');
%       savefig(fullfile(basepath, 'Bead_tracking',['All-Reject-Legit_fractionVtime_fov' num2str(fovn) '.fig']))
%       print('-f',fullfile(basepath, 'Bead_tracking',['All-Reject-Legit_fractionVtime_fov' num2str(fovn) '.png']),'-dpng');


