function [w_r_col,zeta_r_col]=CE(h_receptance_col,D_t,N_hat)

if length(h_receptance_col)<(2*N_hat)^2+1,error('length(h_receptance_col) must be >=(2*N_hat)^2+1'),end
L=floor(length(h_receptance_col)/(2*N_hat));

beta_col=lsqminnorm(hankel(h_receptance_col(1:L),h_receptance_col(L+(1:2*N_hat)-1)),-h_receptance_col(2*N_hat+(1:L)));
z_col=roots(flipud([beta_col(1:2*N_hat);1]));
s_vec=log(z_col)/D_t;

%Sort eigenvalues
[~,Index]=sort(abs(imag(s_vec)));
s_vec=s_vec(Index);

[w_r_col,zeta_r_col]=MDOF_Modal_Param_Visc(s_vec);