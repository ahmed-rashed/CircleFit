function visualizeLocalReceptance(f_local,Receptance_local,circ_prop,ii,jj,ax)

Df=f_local(2)-f_local(1);
if nargin<6
    ax=cla;
end

plot_FRF_Nyq(Receptance_local,[],"\alpha_{"+ii+','+jj+'}',false,'.-k','MarkerSize',10);
title(ax,"$\qquad\qquad\qquad"+min(f_local)+'\leq f\leq'+max(f_local)+'\;,:\Delta f='+Df+'$','interpreter','latex')

for tt=1:1/Df:length(f_local)
    vec=(Receptance_local(tt)-(circ_prop.x_center+1i*circ_prop.y_center))/20;
    aang=angle(vec);
    axes(ax);
    h=text(real(Receptance_local(tt)-vec),imag(Receptance_local(tt)-vec),"$f="+f_local(tt)+'$','interpreter','latex','FontSize',8);
    if ~(aang>-pi/2 && aang<=pi/2)
        aang=aang-pi;
    else
        set(h,'HorizontalAlignment','right')
    end
    set(h,'Rotation',aang/pi*180)
end