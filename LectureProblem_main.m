clearvars
close all
clc

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1;0,0,0])
set(groot,'DefaultLineLineWidth',1);

% Model parameters:
% mass matrix
M=[1 0 0
    0 1 0
    0 0 1] ;
% damping matrix
C=[50 -20 -20
    -20 50 -30
    -20 -30 50]/7;
% stiffness matrix
K=[237315 -161000 0
    -161000 398315 -161000
    0 -161000 398315];

% FRF storage
f_max=200;
N=400;
D_f=f_max/N;
f_col=(0:N-1).'*D_f;

m_row=[2,3];
n_row=[3,3];

[EigVectors_Normalized,EigValues_vec]=MDOF_Eig_Visc(M,C,K);
Receptance=MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,2*pi*f_col,m_row,n_row);

f_mode_min=[120];
f_mode_max=[125];
ShowInternalDetails=true;
CircleFitProblem(Receptance,D_f,m_row,n_row,f_mode_min,f_mode_max,ShowInternalDetails);

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove')

export_figure(3,'',"CircleFitProblem")