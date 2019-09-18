function renameFOVBF(fovn,filext)

%%  Script to rename image files and make a time vector
%   Name of script is RenameFOV.m
%   the goal is to rename image files with fov prefix that matches the
%   folder name they reside in. Script works on tiff files saves renamed
%   and converted images in a subdirectory
%   Drew Drabek, Blacklow and Loparo Labs, Harvard Medical School, BCMP
%   2016-11-29
%
%   Make sure that there is only one bioformats file in the folder of
%   interest
%
%   Need to have in the MATLAB path the bioformats package
%   This script is based on this webpage: 
%   https://www.openmicroscopy.org/site/support/bio-formats5.5/developers/matlab-dev.html#displaying-images
%
%   2018/2/13 AD - added '-v7.3' to the saving of data to support large file
%   sizes
%% Ensure you are in the working directory with the image of interest
    % Document the working directory for the images
    basepath = [cd,'/'];
    
    % Document the index of the experiment, make subfolder
    % commented out to be used as input for function
    % fovn = 0;
    subfolder = ['fov',num2str(fovn)];
    mkdir(subfolder);
    
    % File extension of raw images
    % commented out to be used as input for function
    % filext = 'tiff';
    
    % find all the file types in the working directory with extension input
    filetype = strcat('*.',filext);
    imagelist = dir(fullfile(filetype));

     % sort imagelist based on natural numbers
    [~,index] = sort_nat({imagelist.name}.');
    imagelist = imagelist(index); 
    clear index;
    
    % Convert the set of image names to a cell array.
    imagenames = {imagelist.name}';
    
    % Read the file into 
    tic
    data = bfopen(imagenames{1});
    toc
    
    % save the parameters and date, version -7.3 is no compression and
    % supports file sizes >=2GB (i.e. most cutting movies)
    save(fullfile(basepath,['fov' num2str(fovn) '_data']),'data','-v7.3');
    
%% Make the fov#_times.mat vector
    % create a filename for the index
    fovtimesname = ['fov',num2str(fovn),'_times.mat'];
    time = (0:1:(length(data{1,1})-1))';
    save(fovtimesname,'time');
    
%%  Start the file copy and rename loop
    % Read the images, set new name, write to output folder
tic
    for i = 1 : length(data{1,1})
        % read the ith image in the list as a 8-bit unsigned
        img = data{1,1}{i,1};
        % index = sscanf(imagenames{i},'%d'); %read  filename for number
        % %comment out from old version from QImage software
        index = time(i);
        
        % create a name for re-saving the file per fov# convention
        output_name = [subfolder,'_',sprintf('%04d',index),'.tif'];
        
        %write the image with the new name into the subfolder specified
        imwrite(img,fullfile(subfolder,output_name));
        
    end
   
toc

%% end