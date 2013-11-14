[data,names] = loadImageData('hw3_dataset/1_2d_3k','png');
options.overlay=0;
data=double(data);
D=L2_distance(data,data,1);
[Y,R,E]=Isomap2(D,data,'k',7,options);

figure;
hold on;
plot(Y.coords{2}(1,:), Y.coords{2}(2,:), 'yo') 
angles = importdata('angles.txt');

A=[7 70 263 524 705 836 913 1273 1593 1816 2332 2507 2816];
for q = 1:size(Y.coords{2},2)
    if ~isempty(find(A==q))
        plot(Y.coords{2}(1,q), Y.coords{2}(2,q), 'bo') 
        str = ['(', num2str(angles(q,1)), ', ', num2str(angles(q,2)), ')'];
        text(Y.coords{2}(1,q),Y.coords{2}(2,q),str, 'Color', 'k','FontSize', 10);
    end
end
title('Variation of Theta1');

figure;
hold on
plot(Y.coords{2}(1,:), Y.coords{2}(2,:), 'yo') 

A=[7 67 129 254 362 697 893 1095 1121 1236 1420 1728 1882 1910 1969 2222 2229 2401 2744];
for q = 1:size(Y.coords{2},2)
    if ~isempty(find(A==q))
        plot(Y.coords{2}(1,q), Y.coords{2}(2,q), 'bo') 
        str = ['(', num2str(angles(q,1)), ', ', num2str(angles(q,2)), ')'];
        text(Y.coords{2}(1,q),Y.coords{2}(2,q),str, 'Color', 'k','FontSize', 10);
    end
end
title('Variation of Theta2');

return;