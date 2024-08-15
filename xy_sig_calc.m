function [xy_sig] = xy_sig_calc(dx, x, dy, y, z)
dx = nonzeros(dx);
x = nonzeros(x);
dy = nonzeros(dy);
y = nonzeros(y);
z = nonzeros(z);

if dx ==0 | x==0 | dy==0 | y==0 | z==0
    xy_sig = 0;
else
    xy_sig = sqrt((dx./x).^2+(dy./y).^2).*z;
    xy_sig = xy_sig';
end