function [vol] = calc_vol(v, gauss_func)
[X,Y] = meshgrid(v,v);
A = gauss_func.a1;
x0 = gauss_func.b1;
y0 = x0;
sigmaX = gauss_func.c1;
sigmaY = sigmaX;
Z = A*exp(-((X-x0).^2./(2*sigmaX^2)+(Y-y0).^2./(2*sigmaY^2)));
C = Z;
% figure(1);
% surf(X,Y,Z,C);  crameri('batlowW');
% colorbar
fun = @(x,y) A*exp(-((x-x0).^2./(2*sigmaX^2)+(y-y0).^2./(2*sigmaY^2)));
vol = integral2(fun,v(1),v(end),v(1),v(end));
end
