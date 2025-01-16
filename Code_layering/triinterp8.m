function [Y] = triinterp8(XYZ,y),

% Linear interpolation of current location
xBas = [0 0 0 0 1 1 1 1] + floor(XYZ(1,1));
yBas = [0 0 1 1 0 0 1 1] + floor(XYZ(1,2));
zBas = [0 1 0 1 0 1 0 1] + floor(XYZ(1,3));
xCom = XYZ(1,1) - floor(XYZ(1,1));
yCom = XYZ(1,2) - floor(XYZ(1,2));
zCom = XYZ(1,3) - floor(XYZ(1,3));

% Linear interpolation percentages.
perc=[(1-xCom) * (1-yCom) * (1-zCom); (1-xCom) * (1-yCom) * zCom;
      (1-xCom) *    yCom  * (1-zCom); (1-xCom) *    yCom  * zCom;
         xCom  * (1-yCom) * (1-zCom);    xCom  * (1-yCom) * zCom;
         xCom  *    yCom  * (1-zCom);    xCom  *    yCom  * zCom;];

% The gradient/fiber direction is determined from a voxel
% neighborhood of 8.
Y = 0;
try
    for i=1:8,
        Y  = Y + squeeze(y(xBas(i),yBas(i),zBas(i),:))*perc(i);
    end
catch
    Y = zeros(1,size(y,4));
end
Y = Y(:)';