function [w_r_col,zeta_r_col,stabilized_f_zeta_r_col,A_r_cols]=CE(h_receptance_col,f_col,N_hat_max, ...
                                                                    f_r_tol,zeta_r_tol)   %Optional arguments

if ~iscolumn(f_col),error('f_col must be a column vector.'),end
N_f_max=length(f_col);
if any(size(H_Receptance_col)~=size(f_col)),error('Size of H_Receptance_col must be the same as f_col.'),end

if nargin<4 % Natural frequency tolerance
    f_r_tol=1/100;
else
    if isempty(f_r_tol)
        f_r_tol=1/100;
    end
end

if nargin<5 % zeta tolerance
    zeta_r_tol=5/100;
else
    if isempty(zeta_r_tol)
        zeta_r_tol=5/100;
    end
end

% Sampling parameters
f_max=f_col(end);
N_t=2*N_f_max-2;  %This function assumes N_t=even. This is a valid assumption since N_t is usually power of 2.
% N_t=2*N_f_max-1-ceil(rem(N_f_max/2,1));
D_f=f_col(2)-f_col(1);
[D_t,~,~]=samplingParameters_T_N(1/D_f,N_t);

if N_t<(2*N_hat_max)^2+1,error('N_t must be >=(2*N_hat)^2+1'),end
L=floor(length(h_receptance_col)/(2*N_hat_max));

h_receptance_col=ifft(H_Receptance_col,N_t,1,'symmetric')*N_t/2; % impulse response functions
f_r_cell=cell(N_modes_max,1);
zeta_r_cell=f_r_cell;
stabilized_f_r_cell=f_r_cell;
stabilized_f_zeta_r_cell=f_r_cell;
z_r_cell=f_r_cell;
ind_cell=f_r_cell;
InvConditionNumber_col=zeros(N_modes_max,1);
leastSquareError_col=InvConditionNumber_col;
for N_hat=1:N_hat_max
    h_mat=hankel(h_receptance_col(1:L),h_receptance_col(L+(1:2*N_hat)-1));
    h_dash_col=h_receptance_col(2*N_hat+(1:L));
    beta_col=lsqminnorm(h_mat,-h_dash_col);

    % Error calculation 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    InvConditionNumber=1/cond(h_mat);
    
    %    SingVals=svd(h_T_transpose); %Singular normalized  values
    %    err(N_hat)=1/(max(SingVals)/min(SingVals));
    %		- of least squares
    
    leastSquareError=norm(-h_dash_col-h_mat*beta_col);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    z_r_col=roots(flipud([beta_col(1:2*N_hat);1]));
    s_r_col=log(z_r_col)/D_t;
    
    %Sort eigenvalues
    [~,Index]=sort(abs(imag(s_r_col)));
    s_r_col=s_r_col(Index);
    
    % Modal parameters extraction
    [w_r_current_col,zeta_r_current_col]=pole2modal_Visc_cleaned(s_r_col,D_f,f_max);
    ind_cell{N_hat}=ind1_col(ind2_col);

    f_r_current_col=w_r_current_col/2/pi;
    if N_hat==1
        f_r_cell{N_hat}=f_r_current_col;
        zeta_r_cell{N_hat}=zeta_r_current_col;
        stabilized_f_r_cell{N_hat}=false(size(zeta_r_current_col));
        stabilized_f_zeta_r_cell{N_hat}=false(size(zeta_r_current_col));
    else
        [f_r_cell{N_hat},zeta_r_cell{N_hat},stabilized_f_r_cell{N_hat},stabilized_f_zeta_r_cell{N_hat}]=stabilizeModalParameters(f_r_current_col,zeta_r_current_col,f_r_cell{N_hat-1},zeta_r_cell{N_hat-1},f_r_tol,zeta_r_tol);
    end
end