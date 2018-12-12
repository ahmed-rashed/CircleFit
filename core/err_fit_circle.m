function [x_0,y_0,R_0]=err_fit_circle(x,y)

%Based on equation (4.158) of Maia, Silva et. al., "Theoretical and Experimental Modal Analysis", Research Studies Press, 1997.

L=length(x);
SXX=sum(x.^2);
SYY=sum(y.^2);
SXY=sum(x.*y);
SX=sum(x);
SY=sum(y);
SXXX=sum(x.^3);
SYYY=sum(y.^3);
SXXY=sum(x.^2.*y);
SXYY=sum(x.*y.^2);
A=[SXX SXY SX;SXY SYY SY;SX SY L];
B=[-(SXXX+SXYY);-(SXXY+SYYY);-(SXX+SYY)];
sol=A\B;
a=sol(1);
b=sol(2);
c=sol(3);
x_0=-a/2;
y_0=-b/2;
R_0=sqrt(x_0^2+y_0^2-c);
