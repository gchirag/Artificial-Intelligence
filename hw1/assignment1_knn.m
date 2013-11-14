images = loadMNISTImages('train-images.idx3-ubyte');
labels = loadMNISTLabels('train-labels.idx1-ubyte');
sample = loadMNISTImages('t10k-images.idx3-ubyte');
sample_labels = loadMNISTLabels('t10k-labels.idx1-ubyte');
smpl=sample';
img=images';

%count=0;i=1;k=1;j=1;
%c=cell(28,28);
%for i=1:784
%       c(floor((i-1)/28)+1,mod(i-1,28)+1)=img(94,i);
%end
%labels(94)
%imshow(c')
%var1=knnclassify(smpl,img,labels,1);
%ans1=sum(var1-sample_labels(:)==0);
%var2=knnclassify(smpl,img,labels,2);
%ans2=sum(var2-sample_labels(:)==0);
%var3=knnclassify(smpl,img,labels,3);
%ans3=sum(var3-sample_labels(:)==0);
%var4=knnclassify(smpl,img,labels,4);
%ans4=sum(var4-sample_labels(:)==0);
%var5=knnclassify(smpl,img,labels,5);
%ans5=sum(var5-sample_labels(:)==0);
%var6=knnclassify(smpl,img,labels,6);
%ans6=sum(var6-sample_labels(:)==0);
%var7=knnclassify(smpl,img,labels,7);
%ans7=sum(var7-sample_labels(:)==0);
%[images2, labels2] = loadDigits(3000, 'train','./');
%img2=images2';

success_value=ones(1,50);
error_percentage=ones(1,50);
for i=1:50
   var=knnclassify(smpl,img,labels,i);
   success_value(i)=sum(var-sample_labels(:)==0); 
   error_percentage(i)=10000-success_value(i);
   error_percentage(i)=error_percentage(i)/100;
end

plot(error_percentage)
saveas(figure,'Error_in_knn.png');
set(gca,'XTick',1:2:50);
xlabel('K');
ylabel('Error Percentage');
title('Plot of Error percentages using KNN Classifier');
set(findobj(gca,'Type','Line','Color',[0 0 1]),'Color',[1,0,0]);