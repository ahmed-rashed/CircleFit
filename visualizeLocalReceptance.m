function visualizeLocalReceptance(f_local,Receptance_local,circ_prop,ii,jj,ax)

Df=f_local(2)-f_local(1);
if nargin<6
    ax=cla;
end

plot_FRF_Nyq(Receptance_local,[],['\alpha_{',int2str(ii),',',int2str(jj),'}'],false,'.-k','MarkerSize',10);
title(ax,['$\qquad\qquad\qquad',num2str(min(f_local)),'\leq f\leq',num2str(max(f_local)),'\;,:\Delta f=',num2str(Df),'$'],'interpreter', 'latex')

for tt=1:1/Df:length(f_local)
    vec=(Receptance_local(tt)-(circ_prop.x_center+1i*circ_prop.y_center))/20;
    aang=angle(vec);
    axes(ax);
    h=text(real(Receptance_local(tt)-vec),imag(Receptance_local(tt)-vec),['$f=',num2str(f_local(tt)),'$'],'interpreter', 'latex', 'FontSize',8);
    if ~(aang>-pi/2 && aang<=pi/2)
        aang=aang-pi;
    else
        set(h,'HorizontalAlignment','right')
    end
    set(h,'Rotation',aang/pi*180)
end
