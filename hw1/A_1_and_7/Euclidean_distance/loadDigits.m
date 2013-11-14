function [data, labels] = loadDigits(N, dataSet, dir)
%%
% Description:
%   loadDigits loads the first 'N' digit samples in to 'data' 
%
%
% Parameters:
%   N       : the number of samples to load.
%   dataSet : the data set to load from.
%           : can be one of the two values: 'train' and 'test'
%   dir     : data directory; default value is './data'
%
% Return:
%   data    : D x N matrix containing the first N data samples from the
%             given data set. D is the dimensionality of the input space.
%   labels  : 1 x N matrix containing the corresponding labels.
%
%
% Notes:
%   1. Make sure the data files have the same names as in http://yann.lecun.com/exdb/mnist/
%   2. Unzip the data files into './data' in the working directory.
%
%
% Example: 
%   To read the first 2000 samples in to 'data'
%   and the corresponding labels into 'labels'
%
%   >> [data, labels] = loadDigits(2000, 'train');
%
%
%%

    if nargin < 2
        dataSet = 'train';
    end
    
    if nargin < 3
        dir = './data';
    end
    
    if ~strcmp(dataSet, 'train') && ~strcmp(dataSet, 'test')
        disp('Invalid value for dataSet. Loading digits from the training set.');
        dataSet = 'train';
    end
    
    if strcmp(dataSet, 'train')
        file1 = [dir '/train-images.idx3-ubyte'];
        file2 = [dir '/train-labels.idx1-ubyte'];
    else
        file1 = [dir '/t10k-images.idx3-ubyte'];
        file2 = [dir '/t10k-labels.idx1-ubyte'];
    end
    
    f1 = fopen(file1, 'r', 'b');
    f2 = fopen(file2, 'r', 'b');
    
    if f1 == -1 || f2 == -1
        error('Cannot open data files for reading. Please check the paths.');
    end
    
    magicNum1 = fread(f1, 1, 'int32');
    magicNum2 = fread(f2, 1, 'int32');
    if magicNum1 ~= 2051 || magicNum2 ~= 2049
        error('The data files are corrupt. Please download from http://yann.lecun.com/exdb/mnist/');
    end
        
    N1 = fread(f1, 1, 'int32');
    N2 = fread(f2, 1, 'int32');
    
    if N > N1 || N > N2
        N = min(N1, N2);
        fprintf('There are only %d samples in the data set. Returning only %d samples.\n', N, N);
    end
    
    nr = fread(f1, 1, 'int32');
    nc = fread(f1, 1, 'int32');
    dim = nr * nc;
    
    data = zeros(dim, N);
    labels = zeros(1, N);
    
    a=1;
    
    while(a<=N)
        labels(a)= fread(f2, 1, 'uchar');
        data(1:dim, a)= fread(f1, dim, 'uchar');
        if labels(a)==1 || labels(a)==7
           a=a+1;
        end
    end

end