clear all; 
%preparing workspace
close all;
workspace; % confirm that workspace panel is showing
fontSize = 14;
%frac = 
scale = 1.238736405e-6; % um/pixels 
%scale = (100/108)*1e-6; % um/pixel scale = 1.238736405e-6; % um/pixels 
K_B = 1.38064852e-23; %J/K Boltzmann's constant
T = 298.15; %K Temperature = 25 degree centigrade


%user inputs
steps = 2880; 
nhood = 3; %averaging neighborhood
vstart = 1;
vstop = 26; %start and stop times for video
erode = 4; %number of pixels in erosion element

[File_Name,Path_Name] = uigetfile('*.*','Select the video');
f = fullfile(Path_Name,File_Name);% building full pathname

mov = VideoReader(f); %mov is an object

%Determine the height and width of the frames
vidHeight = mov.Height;
vidWidth = mov.Width;
x = double(0); %placeholder when creating structure to hold doubles
xcap = (vidHeight - rem( vidHeight,nhood))/nhood;       
ycap = (vidWidth - rem( vidWidth,nhood))/nhood;  %used in averaging, new image size = ycap-by-xcap

%Create a MATLAB movie structure array, s.
s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),'colormap',[]);

f = steps/2+1; %holds fft array. if steps = 360, f = 181 bins
frames = struct('x',zeros(1,steps,'uint16'),'y',zeros(1,steps,'uint16'),'length',x,'radius',x,'spectrum',zeros(1,f,'double'),'R_avg',x, 'Area',x, 'undulations',zeros (1,steps,'double')); %contains data for every frame

%Read one frame at a time using readFrame until the end of the file is reached. 
%Append data from each video frame to the structure array
k = 1;
mov.CurrentTime = vstart;  %start point

while ( (hasFrame(mov)) && (mov.CurrentTime < vstop)) %assumes that it reads frames in ascending order    
        
        s(k).cdata = readFrame(mov);
        %converting to 2-D grayscale
        s(k).cdata = rgb2gray(s(k).cdata); %**will not work if reading 2D video 
        k = k+1;    
end
k = k - 1; %number of actual frames captured

I = s(1).cdata;


%Getting user input for the approximate center and minimum inner and max 
%outer radius helps with speeding up the boundary extraction and prevents
% running out of bounds

imshow(I); %display first frame in the video
[x,y] = ginput(3);  %room for improvement, to make it user friendly?

x_center = y(1);
y_center = x(1);
x_inner  = y(2);
y_inner  = x(2);
x_outer  = y(3);
y_outer  = x(3);

R_in = sqrt ((x_inner - x_center)^2 + (y_inner - y_center)^2);
R_out = sqrt ((x_outer - x_center)^2 + (y_outer - y_center)^2); 

mask = uint8 (circularMask(vidWidth,vidHeight,x_center,y_center,R_in,R_out));
mask2 = logical(mask);
f = @(x) imadjust(x);
blackWhite = @(x) im2bw(x,0.15);
h = @(x) adapthisteq(x, 'NumTiles',[2 2], 'ClipLimit',0.3);
compl = @(x) imcomplement(x);
blob = @(x) bwareaopen (x,8000);
blank = @(x) x+255;
%boundaries set and center selected
for i = 1:k
    
    %Creating mask to get ROI
   % mask = uint8 (circularMask(vidWidth,vidHeight,x_center,y_center,R_in,R_out));
    I1 =  adapthisteq(s(i).cdata, 'NumTiles',[2 2], 'ClipLimit',0.3);
    se = strel ('disk',erode); %structuring element for erosion
    I = imerode(I1,se);
    mashed = mask.*I;  %intermediary image
  %  mask2 = logical(mask);

    %Defining the functions for ROI filtering below
%     f = @(x) imadjust(x);
%     blackWhite = @(x) im2bw(x,0.15);
%     h = @(x) adapthisteq(x, 'NumTiles',[2 2], 'ClipLimit',0.3);
%     compl = @(x) imcomplement(x);
%     blob = @(x) bwareaopen (x,5000);
%     blank = @(x) x+255;


    %Performing the image analysis sequence 
    %might need to adjust the functions above, esp BW and h, if having
    %contrast issues
    filtered = roifilt2(mashed,mask2,blackWhite);
    filtered = roifilt2(filtered,mask2,compl);
    filteredmask = roifilt2(filtered,mask2,blob);
    removemask = imcomplement(filteredmask);
    finalImage = uint8(filteredmask).*s(i).cdata;
    Im = roifilt2(s(i).cdata,removemask,blank); %final image    
    I = Im;
    %info :x,y coordinates of contour pixels,the total length of the outline and the R
    %rotating image 360 times and getting the grayscale value for each theta
    track = 0;    
    for deg =0:360/steps:360
        %traces outline by getting min gray value
        [tracein] = graytracevid(x_center,y_center,deg,I,R_in,R_out); % call trace 
        track = track+1;
        frames.radius(track) = tracein(1);%*scale;
       %intensity(track)    = tracein(2);
        frames.x(track)      = tracein(3);
        frames.y(track)      = tracein(4);
    end
    
    frames.radius = frames.radius*scale;
    
    %adding up the length of the outline
    len = 0;
    
    for n=1:steps
        
        %length = length + sqrt( double ( (frames.x(mod(n,360)) - frames.x( mod(n+1,360)))^2  +   (frames.y(mod(n,360)) - frames.y( mod((n+1),360)))^2 ));
         len = len + sqrt( double ( (frames.x(n) - frames.x(n+1))^2  +   (frames.y(n) - frames.y(n+1) )^2 ));  
    end
    frames.length = len*scale;    
    
    %_____FILTER RADII OUTLINE HERE_________ : the
  
    
    windowSize = 20;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    x = frames.radius;
    y = filter(b,a,x); %filtered
    N = length(x);
    half = floor(windowSize/2);
    
    %now correcting for the first few values: 1-windowSize :%
    %they aren't well represented in the filter  
    for j = 1:windowSize
        if half>(j-1) %means we've to modulo backwards

            b = half - (j-1);
            a = half+j;
            val_a = mean (x(1,1:b));
            val_b = mean (x(1,(N-(b-1)):N));

            y(j) = (val_a+val_b)/2;
        else %if half <= (j-1)

            a = j-half;
            b = j+half;
            y(j) = mean(x(a:b));
        end
    end
    frames.radius = y;
    
   
    
    frames.Ravg = mean(frames.radius); %mean radius    
    results(i) = frames;    
end %end for i = 1:k


%getting the FFT of each frame relative to the mean average radius over the
%whole frames
Ravg = mean(cat(1,results.Ravg)); % mean radius across whole video, for all frames
T = 2*pi*Ravg;
N = length(results(1).radius);
t = [0:N-1]/N; % defining time
t = t*T; % length at each point
freq = [0:N/2-1]/T; %/ wavenumber

for i=1:k
    frames.undulations = results(i).radius - Ravg;
    
    %filter undulations here
    windowSize = 20;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    x = frames.undulations;
    y = filter(b,a,x); %filtered
    N = length(x);
    half = floor(windowSize/2);
    
    %now correcting for the first few values: 1-windowSize :%
    %they aren't well represented in the filter  
    for j = 1:windowSize
        if half>(j-1) %means we've to modulo backwards

            b = half - (j-1);
            a = half+j;
            val_a = mean (x(1,1:b));
            val_b = mean (x(1,(N-(b-1)):N));

            y(j) = (val_a+val_b)/2;
        else %if half <= (j-1)

            a = j-half;
            b = j+half;
            y(j) = mean(x(a:b));
        end
    end
    undulations = y;  % movement from mean position
    
    results(i).undulations = y;
    
    
   % undulations = results(i).radius - Ravg; % movement from mean position
    A = abs (fft(undulations))/(N/2);      
    A2 = A(1:N/2).^2; %power of the +ve half    
    P_spectra(i,:) = A2;    
end
    


%average of FFT
undulationSpectra = sum(cat(1,P_spectra)) / k ; %concatenating the spectrum array and averaging
figure,loglog(freq,undulationSpectra,'-*','linewidth',1);
ylabel('<|U_q|>^2 (m^{2})');
xlabel('q (m^{-1})');
title('power spectrum of undulations vs q');

power = undulationSpectra .* freq.^3;
figure,loglog((freq),(power),'-*');
%figure,plot(log10(freq),log10(power),'-*');
ylabel('<|U_q|>^2.q^3 (m^{-1})');
xlabel('q (m^{-1})');
title(' <|U_q|>^2.q^3 vs  q');



ds = (sum(cat(1,results.length)) / k )/360; % = (avg L)/360 == L/2pi
f_s = 1/ds;




figure,imshow(s(1).cdata); hold on
plot(results(1).y,results(1).x);
axis equal;

X = 0:1:steps;
figure, plot (X,results(1).radius);
xlabel('\theta');
ylabel('Fluctuations h (\theta)(um)');
str = sprintf('Fluctuations with respect to Angle theta');
title (str);
grid on;

mem_len = 2*pi*Ravg; %L

slope = 0.07; %enter slope here
slope2 = -1.65; %slope for straight fit

modulus = K_B * T / (slope * 4 * mem_len);
modulus_2 = K_B * T/(10^slope2*4 * mem_len);

hold on;
plot(X,results(1).undulations+Ravg);
