%images = loadMNISTImages('train-images.idx3-ubyte');
%labels = loadMNISTLabels('train-labels.idx1-ubyte');
%sample = loadMNISTImages('t10k-images.idx3-ubyte');
%sample_labels = loadMNISTLabels('t10k-labels.idx1-ubyte');
%smpl=sample';
%img=images';

%count=0;i=1;k=1;j=1;
%c=cell(28,28);
%for i=1:784
%       c(floor((i-1)/28)+1,mod(i-1,28)+1)=img(94,i);
%end
%labels(94)
%imshow(c')



[images2, labels2] = loadDigits(3000, 'train','./../../');
%img2=images2';
distance=L2_distance(images2,images2,1);
[Y,R,E,M,N]=Isomap(images2,'abcdefghijklmnorstuwxyz12357!',labels2,distance,'k',3);
