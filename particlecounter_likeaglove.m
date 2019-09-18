%% nonlinear least squares fitting double exponential
% This script will fit bead survival curves to a double exponential with
% nonlinear least squares fitting

%%
function [ fitresult1 , gof1, ci1] = particlecounter_likeaglove(fovn,ver,t1,t2,fr)
%% initialize script with data and folders of interest
% load the data
    basepath = [cd ''];
    load(fullfile(basepath, ['Track_analysis_v' num2str(ver)], ...
        ['Track_analysis_remainingsummary_fov' num2str(fovn) '_v' num2str(ver) '.mat']));
    basepath = [cd ''];

    % make a folder for outputs
mkdir(['Track_fits_fov' num2str(fovn) '_v' num2str(ver)]);

%% name the fit
    name=['NBID_' num2str(fovn) '_v' num2str(ver) '_exp2'];
%% rename data for the inputs to the fitting
    data = remainingsummary{5,1};

%% format the data for analysis
%fr=2;               %frame rate
%t1=1200;            %start frame for fit
%t2=5490;            %end frame for fit
    x=data(:,2)./fr;     %frame number
    yf=data(:,3);        %fraction of beads
    yb=data(:,1);        %number of beads

% specifically load the time ranges for the x and y data 
% ydata1 (fraction) and ydata2 (count) and xdata (seconds) variables are used
    ydata1=yf(t1:t2); 
    ydata2=yb(t1:t2);
    xdata=x(t1:t2);

%% perform the fitting for fraction of beads

% Set up fittype and options.
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' ,...
                'Upper', [Inf,0,Inf,0] );
    opts.Display = 'Off';
    opts.StartPoint = [1 -0.001 1 -1e-05];

% Fit model to data.
    [fitresult1, gof1] = fit( xdata, ydata1, ft, opts )
    ci1 = confint(fitresult1,0.95)    
% Plot fit with data.
    figure( 'Name', name );
    h = plot( fitresult1, x, yf );
    legend( h, 'fraction beads vs. frame number', 'double exponential', 'Location', 'NorthEast' );
    % Label axes
    ylim([0 1])
    xlabel 'time, seconds'
    ylabel 'fraction highly mobile beads remaining'
    grid on

    % save the figures
     savefig(fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp2_frac_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp2_frac_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);
%% perform the fitting for fraction of beads

% Set up fittype and options.
    ft = fittype( 'exp2' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
                'Upper', [Inf,0,Inf,0] );
    opts.Display = 'Off';
    opts.StartPoint = [1 -0.001 1 -1e-05];
    
% Fit model to data.
    [fitresult2, gof2] = fit( xdata, ydata2, ft, opts );
    ci2 = confint(fitresult2,0.95) ;   

% % set up a single exponential fit
%     ft = fittype( 'exp1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.StartPoint = [1 -0.01];
% % Fit model to data.
%     [fitresult4, gof4] = fit( xdata, ydata2, ft, opts );
    
% Plot fit with data.
    figure( 'Name', name );
    hold on
    hh = plot( fitresult2, x, yb );
%     hhh = plot( fitresult4, x, yb );
    legend( hh, 'number of beads vs. time', 'double exponential', 'Location', 'NorthEast' );
    % Label axes
    ylim([0 max(yb)])
    xlabel 'time, seconds'
    ylabel 'number highly mobile beads remaining'
    grid on
    hold off
    % save the figures
    % save the figures
     savefig(fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp2_count_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp2_count_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);
% %% find the x value where the second exponential starts
%     % endexp1 has the yvalue where the second exponential "starts"
%     endexp1=find(round(fitresult1(x),3) == round(fitresult1.c,3));

%% perform the fitting for fraction of beads single exponential

% Set up fittype and options.
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' ,...
                'Upper', [Inf,0] );
    opts.Display = 'Off';
    opts.StartPoint = [1 -0.001];

% Fit model to data.
    [fitresult3, gof3] = fit( xdata, ydata1, ft, opts );
    ci3 = confint(fitresult3,0.95);    

% Plot fit with data.
    figure( 'Name', name );
    h = plot( fitresult3, x, yf );
    legend( h, 'fraction beads vs. frame number', 'single exponential', 'Location', 'NorthEast' );
    % Label axes
    ylim([0 1])
    xlabel 'time, seconds'
    ylabel 'fraction highly mobile beads remaining'
    grid on
%%
    % save the figures
     savefig(fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp1_frac_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp1_frac_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);
%% perform the fitting for fraction of beads single exponential

% Set up fittype and options.
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' ,...
                'Upper', [Inf,0] );
    opts.Display = 'Off';
    opts.StartPoint = [1 -0.0001];
    
% Fit model to data.
    [fitresult4, gof4] = fit( xdata, ydata2, ft, opts );
     ci4 = confint(fitresult4,0.95);   

% Plot fit with data.
    figure( 'Name', name );
    hold on
    hh = plot( fitresult4, x, yb );
    legend( hh, 'number of beads vs. time', 'single exponential', 'Location', 'NorthEast' );
    % Label axes
    ylim([0 max(yb)])
    xlabel 'time, seconds'
    ylabel 'number highly mobile beads remaining'
    grid on
    hold off
    % save the figures
    % save the figures
     savefig(fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp1_count_' name '_fov' num2str(fovn) '.fig']));
     print('-f',fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
         ['TrackFits_exp1_count_' name '_fov' num2str(fovn) '.png']),'-dpng');
            close(gcf);
% %% find the x value where the second exponential starts
%     % endexp1 has the yvalue where the second exponential "starts"
%     endexp1=find(round(fitresult1(x),3) == round(fitresult1.c,3));


%% save the fits
    
save(fullfile(basepath,['Track_fits_fov' num2str(fovn) '_v' num2str(ver)],...
    ['Fits_fov' num2str(fovn) '_ver' num2str(ver) '.mat']),...
    'fitresult1','gof1','fitresult2','gof2','fitresult3','gof3','fitresult4','gof4',...
    'ci1','ci2','ci3','ci4','name','ver','t1','t2','fr','data')
%% performt the fit with a different function
% % create a function for single and double exponential
% exp1=@(a,b)a*exp(b*xdata)-ydata;
% exp2=@(a,b,c,d)(a*exp(b*xdata) + c*exp(d*xdata))-ydata;
% 
% % run nonlinear least squares fitting for exp1 and exp2
% fit1=lsqnonlin(exp1,100)
% fit2=lsqnonlin(exp2,1)
end
%% updates
% 20190117-updated the fit options to set upper bounds for coefficients to
% prohibit decay values above 0