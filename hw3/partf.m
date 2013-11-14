[data,names] = loadImageData('hw3_dataset/1_2d_3k','png');
obj1 = imread('imageObs1.png');

obj1 = imresize(obj1,[100 100]);
pixels=[];
for i = 1:100                                                       %% calculating pixels pertaining to object
	for j = 1:100
		if(obj1(i,j,1)==0 && obj1(i,j,2)==255 && obj1(i,j,3)==0)
			pixels=[pixels [i;j]];
		end
	end
end

	
intersecting=[];
for img = size(data,2):-1:1                                         %% getting indexes of images of arm intersecting with obstacle
    for i = 1:size(pixels,2)
		if(	~(data(pixels(1,i)+(pixels(2,i)-1)*100,img)==0)||~(data(pixels(1,i)+(pixels(2,i)-1)*100 + 10000,img)==0)||~(data(pixels(1,i)+(pixels(2,i)-1)*100 + 20000,img)==0))
            intersecting = [intersecting img];
			break;
		end
	end
end

data = double(data);
options.dims=1:10;
options.overlay=0;

D = L2_distance(data,data,1);
[Y,R,E]= Isomappartf(D ,names, intersecting, data , 'k' , 7 , options);