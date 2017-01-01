function [x0,y0,R0] = err_fit_circle(x,y)

% ------------------   This file is part of EasyMod   ----------------------------
%  Internal function
%
%  Optimization of circle parameters using least squares method.
%
% Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS


L = length(x) ; 
SXX = sum(x.^2) ; 
SYY = sum(y.^2) ; 
SXY = sum(x.*y) ; 
SX = sum(x) ; 
SY = sum(y) ; 
xtemp = x.^2 ; 
ytemp = y.^2 ; 
SXXX = sum(xtemp.*x) ; 
SYYY = sum(ytemp.*y) ; 
SXXY = sum(xtemp.*y) ; 
SXYY = sum(x.*ytemp) ; 
A = [SXX SXY SX ; SXY SYY SY ; SX SY L] ; 
B = [-(SXXX+SXYY) ; -(SXXY+SYYY) ; -(SXX+SYY)] ; 
X = A\B ; 
a = X(1) ; 
b = X(2) ; 
c = X(3) ; 
x0 = -a/2 ; 
y0 = -b/2 ; 
R0 = sqrt(a^2/4+b^2/4-c) ; 

