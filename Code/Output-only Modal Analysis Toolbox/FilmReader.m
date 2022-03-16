function [data, pv, ph] = FilmReader(ReadRange,NoiseVar,ParRead,RootPath)
% Run FilmReader.m 
% (1) select all files that one needs to be read
% (2) wait untill the process is finished and an output 'data' is returned
%
% Obtain the path and designate the filename
if nargin<3
    ParRead = 'True';
    RootPath = 'L:\My Drive\Graduate study\Research\Projects\';
elseif nargin<4
    RootPath = 'L:\My Drive\Graduate study\Research\Projects\';
end
addpath([RootPath,'OS'])
[file, path, selectedfile] = openfiles('*.png*');
FileType = '.png';
FileName = file{1};
FileName = split(FileName,'.');
FileName = FileName{1};
TempFileName = split(FileName,'_'); 
FileName = strjoin(TempFileName(1:end-1), '_');

if nargin == 0
    a=dir([path, '*.png']);
    m = numel(a);
    NoiseVar = 0;
elseif nargin < 2
    m = ReadRange(2);
    NoiseVar = 0;
else
    m = ReadRange(2);
end

tic
% Get the read path
ReadPath = path;
% ReadPath = uigetdir('','Set the path from which the raw figure is read.');
% ReadPath = [ReadPath,'\'];

% Get the save path
% SavePath = uigetdir('','Set the path to which the trimmed figure is saved.');
% SavePath = [SavePath,'\'];

Image_temp = imread(selectedfile{1});
Image_temp = rgb2gray(Image_temp);
[pv, ph] = size(Image_temp);

data = zeros(m, pv*ph);

% ------------------------------------------------------------------------
% TESTING PARALLEL LOADING
% ------------------------------------------------------------------------
if isequal(lower(ParRead),'false')
    for i = 1:m
        FileNumber = i;
        LoadImageName = [FileName,'_',num2str(FileNumber),FileType];
        Image = imread([path,LoadImageName]);
    %     Image = rgb2gray(Image + uint8(255*randn(size(Image))));
        Image = rgb2gray(imnoise(Image,'Gaussian',0, NoiseVar));
        imagesc(Image)
        data(i,:) = reshape(Image, 1, pv*ph);
        pause(0.00001)
        progress_bar(i,m)
    end
elseif isequal(lower(ParRead),'true')
    % Check if the parpool is running; if not, execute:
    if isempty(gcp('nocreate'))
        clusterinfo = parcluster('local');
        Nworkers = clusterinfo.NumWorkers;
        parpool('local', Nworkers)
    end
    parfor i = 1:m
        FileNumber = i;
        LoadImageName = [FileName,'_',num2str(FileNumber),FileType];
        Image = imread([path,LoadImageName]);
        %     Image = rgb2gray(Image + uint8(255*randn(size(Image))));
        Image = rgb2gray(imnoise(Image,'Gaussian',0, NoiseVar));
        imagesc(Image)
        data(i,:) = reshape(Image, 1, pv*ph);
        pause(0.00001)
%         progress_bar(i,m)
    end
end
toc
disp('Films Read Successfully!')