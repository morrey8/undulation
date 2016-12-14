%plots intensity outline from two given points : two images compared

I1 = s(50).cdata; 
I2 = coarsegrainfilter(I1,nhood);

figure, imshow(I1);
sx = ginput(2);
sx=floor(sx);
hold on;
plot(sx(:,1),sx(:,2));
x = sx(1);
y1= sx(2);
y2 = sx(4);
%w = I1( (s(1,1),s(2,1)) : (s(1,2),s(2,2)) );

prof = I1( x,y2:y1);
figure, plot(prof);
title 'Unfltered Image';
xlabel ' Distance along the line profile';
ylabel 'Grayscale Intensity'

prof2 = I2( x,y2:y1);
figure, plot(prof2);
title 'Filtered Image';
xlabel ' Distance along the line profile';
ylabel 'Grayscale Intensity'


% % figure,plot(w);
% 
% prof1 = I1(
% 
% I1(sx(1,1),sx(2,1))
% I1(sx(1,2),sx(2,2))
% for i = 1:(y1-y2)
% d=I1( x,y1+i)
% end