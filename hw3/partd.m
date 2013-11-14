[data,names] = loadImageData('hw3_dataset/1_2d_3k','png');
options.overlay=0;
data=double(data);
D=L2_distance(data,data,1);
[Y,R,E]=Isomap2(D,data,'k',7,options);

figure;
         hold on;
         plot3(Y.coords{3}(1,:), Y.coords{3}(2,:),Y.coords{3}(3,:), 'ro'); 
         box on;
         view(26, 42);
         xlabel('x-axis');
         ylabel('y-axis');
         zlabel('z-axis');

         title('3D plot for 3000 images');

         hold off;
            
             
return;            