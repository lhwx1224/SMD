function [Y, Nx, Ny, fs] = FilmSnap4OMA(Coordinates, SnapRange, Newfps)
% FilmSnap4DIC snaps and resamples a given video file (default is .AVI
% format). This function is created for the purpose of obtaining the ideal
% range and frame rate of a video for digital image correlation.
% Coordinates: there are 3 options - 
%              1. 4*1 array [xRange, yRange], where xRange is a 2*1 array
%                 whose first element is the starting x position and whose
%                 second element is the ending x position, yRange is also a
%                 2*1 array whose first and second entries are the starting
%                 and ending y coordiantes, respectively.
%              2. string - 'pick', manually picking the range of the snap
%                 picking the range manually when the first snap pops up
%                 using the two-point picking method.
%              3. string - 'full', automatically set the range to the full
%                 width and height of the frame.
% SnapRange - 2*1 array [StartNum, EndNum], where the StartNum and EndNum
%             dictate the starting and ending frame #.
% Newfps  -   integer, New frame rate (fps) after the resampling,  that
%             designate the resampling/downsampling rate, e.g., if the
%             oringial sampling rate is 10 fps and one choose Newfps =
%             1, the sampling rate afterwards is set to 1 fps.
%
% Example: say one likes to downsample the frame rate to 1 fps (i.e., 1 Hz)
%          with a full range of snapping and extract only from the first
%          500 frames -> call: FilmSnap4DIC('full', [1 500], 1) -> select
%          the file one likes to snap and resample -> select the folder to
%          which the new files are saved -> if the 'pick' opetion is
%          selected, pick two points that encircles a square area -> the
%          conversion will start automatically.
%
% Hewenxuan Li Aug. 5 2021

% File paths
[file, path] = uigetfile({'*.mp4'; '*.mov'; '*.avi'},'Select the video file','FileName');

SavePath = uigetdir('','Set the path to which the trimmed figure is saved.');
SavePath = [SavePath,'\'];

if isequal(file, 0)
    disp('No file has selected! User canceled the selection!')
else
    selectedfile = fullfile(path,file);
    FilmInfo = VideoReader(selectedfile); % load relative displacement data
    Nframe = FilmInfo.NumFrame;           % Get the total # of frames
    fs = FilmInfo.FrameRate;              % Get the sampling rate
    disp(['Selected File:', selectedfile]) % Display the file name to be loaded.
end
% Designate the range of the plot to crop
if isequal(lower(Coordinates), 'pick')
    clf
    imagesc(read(FilmInfo,1));
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

if nargin < 2
    SnapRange = [1 FilmInfo.NumFrames];
    Newfps =  FilmInfo.FrameRate;
elseif nargin < 3
    Newfps =  FilmInfo.FrameRate;
end

% Time range for snapping
StartNum = SnapRange(1);
EndNum = SnapRange(2);

if EndNum > Nframe
    EndNum = Nframe;
    disp('Requested frame exceeded the maximum frame! Adjusted to the maximum frame.')
end
% Time range setup
Start = StartNum;
End = EndNum;
Resampling = fs/Newfps;
Indx = Start:Resampling:End; % New frame index

FileName = split(file,'.');
FileName = FileName{1};
% FilmWrite = VideoWriter('newmov.avi');
% open(FilmWrite);
Ny = yRange(2) - yRange(1) + 1;
Nx = xRange(2) - xRange(1) + 1;
Y = zeros(Nx*Ny, length(Indx));
% Process the images
for i = 1:length(Indx)
    clf
%     SaveName = [FileName,'_fps',num2str(Newfps),'_image_',num2str(i),'.tiff'];
    Image = read(FilmInfo,i); % Read the resampled frame
    Image_trim = Image(yRange(1):yRange(2), xRange(1):xRange(2), :); % Trim the frame
    Image_trim = rgb2gray(Image_trim);
    imshow(Image_trim)
    colormap('gray')
%     frame = getframe(gcf);
    Y(:,i) = reshape(Image_trim, Nx*Ny, 1);
%     writeVideo(FilmWrite,frame);
%     imwrite(Image_trim, [SavePath, SaveName])
end
% close(FilmWrite);