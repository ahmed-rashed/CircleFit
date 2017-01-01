function CircleFitInternalDetails(f_local,H_local,H_gen_local,circ_prop,theta_horizontal_r,f_before_mat,f_after_mat,eta_r_mat)
figure

% Nyquist's circle visualization
subplot(6,2,[1,3,5]);
set(gca,'XAxisLocation','top')
plot_FRF_Nyq([H_local,H_gen_local],'H^\textrm{local}',false,'.-');

% w_r point on the Nyquist diagram
x_r = circ_prop.x_center+circ_prop.Radius*cos(theta_horizontal_r);
y_r = circ_prop.y_center+circ_prop.Radius*sin(theta_horizontal_r);
hold on
plot([circ_prop.x_center,x_r],[circ_prop.y_center,y_r],'k');

% FRF magnitude
ax_mag_h=subplot(6,2,[2,4,6,8]);
ax_phase_h=subplot(6,2,[10,12]);
plot_FRF_mag_phase(f_local,[H_local,H_gen_local],false,ax_mag_h,ax_phase_h,'','H^\textrm{local}');

legend('Measured FRF','Generated FRF');

subplot(6,2,[7,9,11]);
surf(f_before_mat,f_after_mat,eta_r_mat);
axis tight
xlabel('$f_{\textrm{before}}$, (Hz)','interpreter', 'latex');
ylabel('$f_{\textrm{after}}$, (Hz)','interpreter', 'latex');
zlabel('$\eta_r$','interpreter', 'latex');