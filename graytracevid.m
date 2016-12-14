function [out] = graytracevid (x0,y0,theta,I,Rin,Rout)
% tuesday 10/04/2016
%takes angle and moves outwards checking for minimum value between specified radii
xe = x0;  %input arguments
ye = y0;
%Rin  = round (Rin-0.5);
%Rout = round (Rout+0.5);  
%[ytframe , xrframe] = size(I);

[xrframe,ytframe] = size(I);
tol = 10; % border tolerance
step = Rin;
%brun = 1;

% Preallocating arrays for efficiency. 
if xrframe > ytframe
    intensityProfile = zeros(1,xrframe);
    xaxis   = zeros(1,xrframe);
elseif xrframe < ytframe
    intensityProfile = zeros(1,ytframe);
    xaxis   = zeros(1,ytframe);
else
    intensityProfile = zeros(1,xrframe);
    xaxis   = zeros(1,xrframe);
end

i = 1 ;

xcoord = round(xe+step*cosd(theta));
ycoord = round(ye+step*sind(theta));

while ~(xcoord >= (xrframe - tol) || xcoord <= (0+tol)|| ycoord >= (ytframe - tol) || ycoord <= (0 + tol) || step >= Rout)
    xcoord = round(xe+step*cosd(theta));
    ycoord = round(ye+step*sind(theta));
    intensityProfile(i) = I(xcoord,ycoord);
    xaxis(i) = step; 
    i = i + 1;
    step = step + 1;
end
% shrinking the arrays
intensityProfile = intensityProfile(1:i-1);
xaxis = xaxis (1:i-1);


minIntesity = min(intensityProfile);
% s = size(intensityProfile); % could use step or i-1
%you may also use i-1 instead of step : (update - I actually did!)
ind = 1;
for j = 1 : i-1 %finding the indeces @ which min happens
%checking if more than one index location has min
    if intensityProfile(j) <= minIntesity
        index(ind) = j;
        ind = ind + 1;
    end
end

%{ 
%uncomment to take care of multiple min values
plot(xaxis,intensityProfile);
k = 1;
n = ind -1; %size(index);
sum = [0.0;0.0];
%  may also use ind-1 instead of n
if n > 1 %more than 1 minimums? && remember to add a control to check the standard deviation
%     for k = 1 : n
%         q = index(k);
%         sum (1,k+1) = sum(1,k) + intensityProfile(q); %summing up the intensities about indeces in index
%         sum (2,k+1) = sum(2,k) + index(k); % summing up the indeces in order to  find distance travelled. is it equal to step?
%         k = k +1;
%     end
out(1:4) = 0; %sum/n; % length to outer edge
%out(2) = 0; %sum/n; % average intensity


else

%}
    out(1) = index(1)+Rin; %length to outer edge (step)
    out(2) = minIntesity; %intensity
    out(3) = x0+(index(1)+Rin)*cosd(theta); % x
    out(4) = y0+(index(1)+Rin)*sind(theta); % y
    
end
%end

