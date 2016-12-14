function coarse = coarsegrainfilter (I,nhood) 
% replaces all pixels in nhood by nhood box with average value
s = size(I);
height = s(1);
width  = s(2);
coarse = 255*ones(height,width); % will hold processed image

%user inputs
nhood = 3; % pixels to include in the averaging

h = floor(height/nhood); %number of boxes stacked upwards in new image (coarse)
b = floor(width/nhood); %number of boxes stacked sideways in new image (coarse)
h_rem = mod(height,nhood); % if I need it later for some reason
b_rem = mod(width,nhood); % if I need it later for some reason

 for x = 0:h-1        
        for y = 0:b-1
            
            sum = double(0);
            for n = 1:nhood                
                for m = 1:nhood                    
                    sum = double (sum + double(I(x*nhood+n , y*nhood+m))); %sum of pixels in neighborhood
                end
            end
            coarse(x*nhood+1:(x*nhood+nhood),y*nhood+1:(y*nhood+nhood)) = uint8 (sum / nhood^2);    % in
        end
 
 end
 coarse = uint8(coarse);
%improfile: plots intensity along chosen line
% x = [19 427 416 77];
% y = [96 462 37 33];
% improfile(I,x,y),grid on;
