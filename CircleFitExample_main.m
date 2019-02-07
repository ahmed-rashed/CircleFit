clearvars
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

% Model parameters:
% mass matrix
M=[1 0 0; 0 1 0; 0 0 1] ;
% damping matrix
C=[40 0 0; 0 40 0; 0 0 40] ;
% stiffness matrix
K=[237315 -161000 0; -161000 398315 -161000; 0 -161000 398315] ;

% FRF storage
f_max=200;
N=400;
D_f=f_max/N;
f_col=(0:N-1).'*D_f;

N_DOF=size(M,1);
n_row=ones(1,N_DOF);
m_row=1:N_DOF;

[EigVectors_Normalized, EigValues_vec]=MDOF_Eig_Visc(M, C, K);
Receptance=MDOF_FRF_Visc(EigValues_vec, EigVectors_Normalized, 2*pi*f_col, n_row, m_row);

f_mode_min=[45 85 115] ;
f_mode_max=[55 95 125] ;
ShowInternalDetails=true;
CircleFitProblem(Receptance,D_f,n_row,m_row,f_mode_min,f_mode_max,ShowInternalDetails);

[w_r_col_exact, zeta_r_col_exact]=MDOF_Modal_Param_Visc(EigValues_vec);
w_r_col_exact/2/pi, zeta_r_col_exact

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove')

export_figure((1:3),'==',{'CircleFit1','CircleFit3','CircleFit2'})
NN=N_DOF*length(f_mode_min);
file_names=cell(1,NN);
for ii=1:NN
    file_names{ii}=['InternalDetails',int2str(ii)];
end
export_figure(3+(1:NN),'==',file_names)
