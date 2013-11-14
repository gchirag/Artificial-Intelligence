[data,names] = loadImageData('hw3_dataset/1_2d_3k','png');

intersecting=[];

data = double(data);
options.dims=1:10;
options.overlay=0;

D = L2_distance(data,data,1);
[Y,R,E]= Isomappath(D ,names, intersecting, data , 'k' , 7 , options);