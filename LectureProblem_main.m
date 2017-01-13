clearvars
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

% Model parameters:
% mass matrix
M = [1 0 0
    0 1 0
    0 0 1] ;
% damping matrix
C = [50 -20 -20
    -20 50 -30
    -20 -30 50];
% stiffness matrix
K = [237315 -161000 0
    -161000 398315 -161000
    0 -161000 398315] ;

% FRF storage
f_max=200;
N=400;
D_f=f_max/N;
f_col=(0:N-1).'*D_f;

n_row=[1,1];
m_row=[1,3];

[EigVectors_Normalized, EigValues_mat]=MDOF_Eig_Visc(M, C, K);
Receptance=MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, 2*pi*f_col, n_row, m_row);

f_mode_min=[49];
f_mode_max=[55];
ShowInternalDetails=true;
CircleFitProblem(Receptance,D_f,n_row,m_row,f_mode_min,f_mode_max,ShowInternalDetails);

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove')

export_figure(3,'',{'CircleFitProblem'})