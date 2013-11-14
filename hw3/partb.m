[data,names] = loadImageData('hw3_dataset/1_2d_3k','png');
options.overlay=0;
data=double(data);
D=L2_distance(data,data,1);
[Y,R,E]=Isomap2(D,data,'k',7,options);

figure;
         hold on;
         plot(Y.coords{2}(1,:), Y.coords{2}(2,:), 'ro'); 
         

         for q = 1:size(Y.coords{2},2)
                if mod(q,150)==0 && ~(q==2100) && ~(q==2550) && ~(q==1500) && ~(q==900) && ~(q==750)
                img=reshape(data(:,q),100,100,3);
               
                scalex = range(Y.coords{2})/800;
                scaley = range(Y.coords{2})/800;
                  
                img = img(21:80,21:80,:);

                img=img/255;
                xc=Y.coords{2}(1,q);
                yc=Y.coords{2}(2,q);

                imagesc([xc xc-100*2*scalex],[yc yc-100*2*scaley], img);
                end
         end
         
         title('Two-dimensional Isomap.');
         hold off;
		 
return;