function [xlocation,ylocation]=LED_location(xstart,ystart,arraysize)
% LED lighting sequence spirals out
% x+ y+ x- x- y- y- x+ x+ x+ y+ y+ y+ 
% <<  <<  .  .  .  .   .
% \/  .  16 15 14  13  .
% \/  .   5  4  3  12  . 
% \/  .   6  1  2  11  .
% \/  .   7  8  9  10  .
% \/  .   .  .  .  .   .
% \/
% Input: 
%       xstart,ystart: absolute coordinate of initial LED
%       arraysize: side length of lit LED array
% Output:
%       xlocation,ylocation: absolute coordinate of lit LEDs

xorynode=zeros(1,70); % satisfy the 32*32 LED matrix
xorynode(1,1)=1;
dif=1;
dif_judge=1;
for i=2:70
    xorynode(1,i)=xorynode(1,i-1)+dif;
    if (dif_judge<2)
        dif_judge=dif_judge+1;
    else
        dif=dif+1;
        dif_judge=1;
    end
end
% light LEDs with NO LED in the center
xlocation=zeros(1,arraysize.^2);
ylocation=zeros(1,arraysize.^2);
xlocation(1,1)=xstart;
ylocation(1,1)=ystart;
xy_order=2;
for i=2:arraysize^2
    if(i<=xorynode(1,xy_order))
    else
        xy_order=xy_order+1;
    end
    if ((rem(xy_order,2))==0)
        xlocation(1,i)=xlocation(1,i-1)+((-1).^(rem(xy_order/2,2)+1));
        ylocation(1,i)=ylocation(1,i-1);
    elseif ((rem(xy_order,2))==1)
        xlocation(1,i)=xlocation(1,i-1);
        ylocation(1,i)=ylocation(1,i-1)+((-1).^(rem((xy_order-1)/2,2)+1));
    end
end

end