function [vol, dI] = calc_vol(v, gauss_func, ci, fA, fx0, fy0, fsigmaX, fsigmaY)
[X,Y] = meshgrid(v,v);
A = gauss_func.a1;
A_ci = ci(:,1);
A_std = sqrt(length(v))*(A_ci(2)-A_ci(1))/3.92;
x0 = gauss_func.b1;
x0_ci = ci(:,2);
x0_std = sqrt(length(v))*(x0_ci(2)-x0_ci(1))/3.92;
y0 = x0;
y0_std = x0_std;
sigmaX = gauss_func.c1;
sigmaY = sigmaX;
sigma_ci = ci(:,3);
sigmaX_std = sqrt(length(v))*(sigma_ci(2)-sigma_ci(1))/3.92;
sigmaY_std = sigmaX_std;
Z = A*exp(-((X-x0).^2./(2*sigmaX^2)+(Y-y0).^2./(2*sigmaY^2)));
C = Z;
% figure(1);
% surf(X,Y,Z,C);  crameri('batlowW');
% colorbar
fun = @(x,y) A*exp(-((x-x0).^2./(2*sigmaX^2)+(y-y0).^2./(2*sigmaY^2)));
vol = integral2(fun,v(1),v(end),v(1),v(end));

syms A a b c d x0 y0 sigmaX sigmaY
dfdA = double(subs(fA, {A,a,b,c,d,x0,y0,sigmaX,sigmaY}, ...
    {gauss_func.a1, v(1), v(end), v(1), v(end), ...
    gauss_func.b1, gauss_func.b1, gauss_func.c1, gauss_func.c1}));

dfdx0 = double(subs(fx0, {A,a,b,c,d,x0,y0,sigmaX,sigmaY}, ...
    {gauss_func.a1, v(1), v(end), v(1), v(end), ...
    gauss_func.b1, gauss_func.b1, gauss_func.c1, gauss_func.c1}));

dfdy0 = double(subs(fy0, {A,a,b,c,d,x0,y0,sigmaX,sigmaY}, ...
    {gauss_func.a1, v(1), v(end), v(1), v(end), ...
    gauss_func.b1, gauss_func.b1, gauss_func.c1, gauss_func.c1}));

dfdsigmaX = double(subs(fsigmaX, {A,a,b,c,d,x0,y0,sigmaX,sigmaY}, ...
    {gauss_func.a1, v(1), v(end), v(1), v(end), ...
    gauss_func.b1, gauss_func.b1, gauss_func.c1, gauss_func.c1}));

dfdsigmaY = double(subs(fsigmaY, {A,a,b,c,d,x0,y0,sigmaX,sigmaY}, ...
    {gauss_func.a1, v(1), v(end), v(1), v(end), ...
    gauss_func.b1, gauss_func.b1, gauss_func.c1, gauss_func.c1}));

dI = (sqrt( (dfdA*A_std)^2 + (dfdx0*x0_std)^2 + (dfdy0*y0_std)^2 + (dfdsigmaX*sigmaX_std)^2 + ...
     (dfdsigmaY*sigmaY_std)^2));

end
