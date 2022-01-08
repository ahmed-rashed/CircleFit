clearvars
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

% Model parameters:
% mass matrix
M_mat=[1 0 0; 0 1 0; 0 0 1] ;
% damping matrix
C_mat=[40 0 0; 0 40 0; 0 0 40] ;
% stiffness matrix
K_mat=[237315 -161000 0; -161000 398315 -161000; 0 -161000 398315] ;

% FRF storage
f_max=200;
N_f_max=400;
D_f=f_max/N_f_max;
f_col=(0:N_f_max-1).'*D_f;

K=round(2.58*N_f_max);
f_s=K*D_f;
D_t=1/f_s;
t_col=(0:K-1).'*D_t;


N_DOF=size(M_mat,1);
m_row=ones(1,N_DOF);
n_row=1:N_DOF;

[EigVectors_Normalized,EigValues_vec]=MDOF_Eig_Visc(M_mat,C_mat,K_mat);
[w_r_col_exact,zeta_r_col_exact]=MDOF_Modal_Param_Visc(EigValues_vec);
w_r_col_exact/2/pi,zeta_r_col_exact

H_receptance_cols=MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,2*pi*f_col,m_row,n_row);
h_receptance_cols=MDOF_IRF_Visc(EigValues_vec,EigVectors_Normalized,t_col,m_row,n_row);

N_hat=N_DOF;
for nn=1:N_DOF
    [w_r_col,zeta_r_col]=CE(h_receptance_cols(:,nn),D_t,N_hat);
    error=norm((w_r_col-w_r_col_exact)./w_r_col_exact) %#ok<NOPTS,NASGU> 
    error=norm((zeta_r_col-zeta_r_col_exact)./zeta_r_col_exact) %#ok<NOPTS> 
end

f_mode_min=[45 85 115] ;
f_mode_max=[55 95 125] ;
ShowInternalDetails=true;
CircleFitProblem(H_receptance_cols,D_f,m_row,n_row,f_mode_min,f_mode_max,ShowInternalDetails);

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove')

export_figure((1:3),'==',["CircleFit1","CircleFit3","CircleFit2"])
NN=N_DOF*length(f_mode_min);
file_names=strings(1,NN);
for ii=1:NN
    file_names(ii)="InternalDetails"+ii;
end
export_figure(3+(1:NN),'==',file_names)
