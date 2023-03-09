function [w_r_clean_col,zeta_r_clean_col,ind_clean]=pole2modal_Visc_cleaned(s_r_col,D_f,f_max)

[w_r_raw_col,zeta_r_raw_col]=pole2modal_visc(s_r_col);

ind1_col=find((w_r_raw_col>(5*2*pi*D_f)) & (w_r_raw_col<2*pi*f_max) & (zeta_r_raw_col>0));
[w_r_clean_col,ind2_col]=uniquetol(w_r_raw_col(ind1_col),1e-7/max(abs(w_r_raw_col)));    % Eliminate repeated natural frequencies (resulting from conjugate poles) and sorts them
ind_clean=ind1_col(ind2_col);
zeta_r_clean_col=zeta_r_raw_col(ind1_col(ind2_col));