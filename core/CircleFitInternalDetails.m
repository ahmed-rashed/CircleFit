function CircleFitInternalDetails(f_local_vec,alpha_local_vec,alpha_gen_local_vec,circ_prop,f_r,beta_r,f_before_mat,f_after_mat,eta_r_mat,label_str)
figure

Df=f_local_vec(2)-f_local_vec(1);

% Nyquist's circle visualization
subplot(6,2,[1,3,5]);
set(gca,'XAxisLocation','top')
plot_FRF_Nyq([alpha_local_vec,alpha_gen_local_vec],[],label_str,false,'.-');
title(['$',num2str(min(f_local_vec)),'\leq f\leq',num2str(max(f_local_vec)),'\;,:\Delta f=',num2str(Df),'$'],'interpreter', 'latex')

% w_r point on the Nyquist diagram
x_r=circ_prop.x_center+circ_prop.Radius*cos(beta_r);
y_r=circ_prop.y_center+circ_prop.Radius*sin(beta_r);
hold on
plot([circ_prop.x_center,x_r],[circ_prop.y_center,y_r],'k');

% Bode plot
ax_mag_h=subplot(6,2,[2,4,6,8]);
ax_phase_h=subplot(6,2,[10,12]);
[ax_mag_h,ax_phase_h]=plot_FRF_mag_phase(f_local_vec,[alpha_local_vec,alpha_gen_local_vec],false,ax_mag_h,ax_phase_h,'',label_str);
temp1=sort([get(ax_phase_h,'XTick'),f_r]);
ind_r=find(temp1==f_r);
set(ax_phase_h,'XTick',temp1);
set(ax_mag_h,'XTick',temp1);
temp2=get(ax_phase_h,'XTickLabel');
temp2{ind_r}=['{\itf_{r}}=',num2str(f_r,'%.2f' )];
set(ax_phase_h,'XTickLabel',temp2);
ax_phase_h.XAxis.FontName='Times';

legend(ax_mag_h,'Measured','Circle-Fit','Location','southeast');

if all(size(f_before_mat)>1)
    subplot(6,2,[7,9,11]);
    surf(f_before_mat,f_after_mat,eta_r_mat);
    axis tight
    xlabel('$f_{\mathrm{b},i}$ (Hz)','interpreter', 'latex');
    ylabel('$f_{\mathrm{a},j}$ (Hz)','interpreter', 'latex');
    zlabel('$\eta_{r,ij}$','interpreter', 'latex');
else
    warning('f_before_mat is not a matrix. damping surface was skipped')
end
