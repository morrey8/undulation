function mask = circularMask (width,height,X,Y,R_in,R_out)

%inputs:imageSizeX,imageSizeY,X,Y,R_in,R_out
%Gives a circular mask to blank out pixels away from liposome outline

[rowsInImage,columnsInImage] = meshgrid(1:width, 1:height);
% Next create the circle in the image.
innerCircle = (rowsInImage - Y).^2 + (columnsInImage - X).^2 <= R_in.^2;
outerCircle = (rowsInImage - Y).^2 + (columnsInImage - X).^2 <= R_out.^2;

% innerCircle and OuterCircle are 2D "logical" arrays.
% Now, display it.
% image(innerCircle);
% image(outerCircle) ;
% colormap([0 0 0; 1 1 1]);
% title('Binary image of the circles');
mask = innerCircle + outerCircle ;
%mask = uint8(mask);
n = size (mask);

for i = 1: n(1)
    for j = 1: n(2)
        
        if (mask(i,j) == 1)
            mask(i,j) = 1;
        else
            mask(i,j) = 0;
        end
    end
end
%figure,imshow(mask);
%m = mashed;
%h =  adapthisteq(mask, 'NumTiles',[28 28], 'ClipLimit',0.1); figure,imshow(h)

