function FilmSnapper(Coordinates, SnapRange, Resampling, Parallel, RootPath)
% FilmSnapper snaps/trims the image file from the acquired motion pictures
% either generated from the numerical simulations or from recorded videos.
% syntax: FilmSnapper(Coordinates, SnapRange, Resampling)
% input: Coordinates - (1) input x range and y range for the program to
%        snap the film [xRange, yRange] or (2) 'pick' picking the
%        upper-left and lower-right corners with an example picture.
%        SnapRange - [StartNum, EndNum] time stamp in between which the
%        film will be included in the trimming process
%        Resamping - integer, Resamping rate for down-sampling the film
%        Parallel - 
% Hewenxuan Li, November 2021
% Output-only Modal Analysis Toolbox v0.0
if nargin<5
    RootPath = 'L:\My Drive\Graduate study\Research\Projects\';
    Parallel = 'True';
elseif nargin<4
    Parallel = 'True';
end

addpath([RootPath,'OS'])
disp('Select the image file to locate the path and to extract the file name:')
% Obtain the path and designate the filename
[file, path, selectedfile] = openfiles('*.*');
FileType = '.png';
FileName = file{1};
FileName = split(FileName,'.');
FileName = FileName{1};
TempFileName = split(FileName,'_'); 
FileName = strjoin(TempFileName(1:end-1), '_');
disp('File location and name are extracted!')
% Get the read path
ReadPath = path;
% ReadPath = uigetdir('','Set the path from which the raw figure is read.');
% ReadPath = [ReadPath,'\'];

disp('Select the folder to which the trimmed figures will be saved:')
% Get the save path
SavePath = uigetdir('','Set the path to which the trimmed figure is saved.');
SavePath = [SavePath,'\'];

% xRange = Coordinates(1:2);
% yRange = Coordinates(3:4);
StartNum = SnapRange(1);
EndNum = SnapRange(2);

Start = StartNum;
End = EndNum;
Indx = Start:Resampling:End;

% Designate the range of the plot to crop
if isequal(lower(Coordinates), 'pick')
    clf
    Image = imread(selectedfile{1});
    imagesc(Image)
    grid minor
    [x,y] = ginput(2);
    xRange = ceil(x);
    yRange = ceil(y);
    hold on, plot([xRange(1) xRange(1)],[yRange(1) yRange(2)],'r')
    plot([xRange(2) xRange(2)],[yRange(1) yRange(2)],'r')
    plot([xRange(1) xRange(2)],[yRange(1) yRange(1)],'r')
    p1 = plot([xRange(1) xRange(2)],[yRange(2) yRange(2)],'r');
    legend(p1,'Snapped Range')
    disp(['Selected Region (in pixel):'])
    disp(['x Range: ', num2str(xRange(1)),' - ',num2str(xRange(2))]);
    disp(['y Range: ', num2str(yRange(1)),' - ',num2str(yRange(2))]);
elseif isequal(lower(Coordinates), 'full')
    xRange = [1 FilmInfo.Width];
    yRange = [1 FilmInfo.Height];
else
    xRange = Coordinates(1:2);
    yRange = Coordinates(3:4);
    disp(['Selected Region (in pixel):'])
    disp(['x Range: ', num2str(xRange(1)),' - ',num2str(xRange(2))]);
    disp(['y Range: ', num2str(yRange(1)),' - ',num2str(yRange(2))]);
end

if isequal(lower(Parallel),'false')
    tic
    for i = 1:length(Indx)
        FileNumber = num2str(Indx(i));
        NewFileNumber = num2str(i);
        LoadImageName = [FileName,'_',FileNumber,FileType];
        SaveName = ['trim','_fr_',num2str(Resampling),'_', FileName, '_', NewFileNumber,FileType];
        Image = imread([ReadPath,LoadImageName]);
        Image_trim = Image(yRange(1):yRange(2), xRange(1):xRange(2), :);
        imwrite(Image_trim, [SavePath, SaveName])
        progress_bar(i, length(Indx))
    end
    toc
elseif isequal(lower(Parallel),'true')

    % Check if the parpool is running; if not, execute:
    if isempty(gcp('nocreate'))
        clusterinfo = parcluster('local');
        Nworkers = clusterinfo.NumWorkers;
        parpool('local', Nworkers)
    end
    tic
    parfor i = 1:length(Indx)
        FileNumber = num2str(Indx(i));
        NewFileNumber = num2str(i);
        LoadImageName = [FileName,'_',FileNumber,FileType];
        SaveName = ['trim','_fr_',num2str(Resampling),'_', FileName, '_', NewFileNumber,FileType];
        Image = imread([ReadPath,LoadImageName]);
        Image_trim = Image(yRange(1):yRange(2), xRange(1):xRange(2), :);
        imwrite(Image_trim, [SavePath, SaveName])
%         progress_bar(i, length(Indx))
    end
    toc
end
finishbeep()