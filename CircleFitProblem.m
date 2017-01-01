function CircleFitProblem(Receptance,Df,n_row,m_row,f_mode_min,f_mode_max,ShowInternalDetails)

N=size(Receptance,1);
f_col=(0:N-1).'*Df;
n_FRF=length(n_row);
n_modes=length(f_mode_min);

% FRF Visualization
f0=figure;
f1=figure;
for ii=1:n_FRF
    for fig=[f0,f1]
        figure(fig)
        ax_mag_h=subplot(2,n_FRF,ii);hold on
        
        label_str=['\alpha_{',int2str(n_row(ii)),int2str(m_row(ii)),'}'];
        plot_FRF_mag_phase(f_col,Receptance(:,ii),false,ax_mag_h,gobjects,'',label_str);

        subplot(2,n_FRF,n_FRF+ii);hold on
        %coloured_line_3d(real(Receptance(:,ii)),imag(Receptance(:,ii)),zeros(size(Receptance(:,ii))),f_col);view(2)
        plot_FRF_Nyq(Receptance(:,ii),label_str);
    end
end

% Circle-fit
f2=figure;
f_r_calc_col=nan(n_modes,n_FRF);
zeta_r_calc_col=nan(n_modes,n_FRF);
for ii=1:n_FRF
    Receptance_Calculated=zeros(N,1);
    for jj=1:n_modes
        Receptance_temp=Receptance(:,ii);
        LocalZone_flag=(f_col>=f_mode_min(jj)) & (f_col<=f_mode_max(jj));
        Receptance_local = Receptance_temp(LocalZone_flag);
        freq_local = f_col(LocalZone_flag);
        
        %Circle Fit
        [Mode_Inf,circ_prop]=FRF_CircleFit(freq_local,Receptance_local,ShowInternalDetails)
        f_r_calc_col(jj,ii)=Mode_Inf.f_r;
        zeta_r_calc_col(jj,ii)=Mode_Inf.eta_r;
        
        Receptance_Calculated=Receptance_Calculated+Mode_Inf.A_r./(complex((2*pi*Mode_Inf.f_r)^2-(2*pi*f_col).^2,Mode_Inf.eta_r*(2*pi*Mode_Inf.f_r)^2));
        
        %Receptance_local visualization
        figure(f2)
        ax=subplot(n_FRF,n_modes,(ii-1)*n_modes+jj);
        visualizeLocalReceptance(freq_local,Receptance_local,circ_prop,n_row(ii),m_row(ii),ax);
    end
    
    figure(f1)
    subplot(2,n_FRF,ii)
    semilogy(f_col,abs(Receptance_Calculated))
    legend('Measured','Calculated');
   
    subplot(2,n_FRF,n_FRF+ii)
    %coloured_line_3d(real(Receptance(:,ii)),imag(Receptance(:,ii)),zeros(size(Receptance(:,ii))),f_col)
    plot(real(Receptance_Calculated),imag(Receptance_Calculated))
    legend('Measured','Calculated');
end