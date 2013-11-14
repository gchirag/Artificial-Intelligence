function [ data, names ] = loadImageData( srcDir, type)
%%  Description: loadImageData
%       - loads all the images of the given directory into a matrix
%
%   Parameters:
%       - srcDir: Source directory
%       - type  : File type. (for example jpg)
%               : = jpg if nothing is given
%
%   Output:
%       - data  : Output matrix containing all the images 
%               : as 10,000 dimensional vectors
%       - names : File names corresponding to each image.
%
%
%   Example: To load all the images from 'nao200' directory
%      
%       >> [data, names] = loadImageData('nao200', 'jpg');
%
%
    fprintf('Loading image data');
    if(nargin < 2)
        type = 'jpg';
    end
    data = [];
    names = {};
    files = dir([srcDir '/*.' type]);
    
    for i=1:length(files)
        if(mod(i, 50) == 0)
            fprintf('.');
        end
        img = imread([srcDir '/' files(i).name]);
        img = imresize(img, [100 100]); 
        [irow icol] = size(img);

        temp = reshape(img,irow*icol,1);    % Reshaping 2D images into 1D image vectors
        data = [data temp];                 % 'data' grows after each iteration
        temp = regexp(files(i).name, '\.', 'split');
        names{i} = temp(1);
    end
    fprintf('done!\n');
end