function [f_r,eta_r,A_r,B_r,circ_prop]=FRF_CircleFit(f_local_vec,alpha_local_vec,ShowInternalDetails,label_str)

if nargin<4
    label_str='\alpha';
else
    if isempty(label_str)
        label_str='\alpha';
    end
end

%Error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the case of insufficient number of samples
N_pts=length(f_local_vec); 
if N_pts < 5
    error('Insufficient number of points in the studied frequency range');
end

%Check f_local_vec is equispaced with a positive d_f
d_f=f_local_vec(2)-f_local_vec(1);
n=(f_local_vec(2:end)-f_local_vec(1:end-1))/d_f;
if (any(abs(n-1)>eps*1000))
    error('Input frequency vector must be equispaced with a positive d_f');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Circle Fitting
[x_center,y_center,Radius]=err_fit_circle(real(alpha_local_vec),imag(alpha_local_vec));
circ_prop=struct('x_center',x_center,'y_center',y_center,'Radius',Radius);

%Identify the two points surrounding w_r
chord=sqrt((real(alpha_local_vec(2:end))-real(alpha_local_vec(1:end-1))).^2+(imag(alpha_local_vec(2:end))-imag(alpha_local_vec(1:end-1))).^2);
[chord_max,index]=max(chord);
f_b=f_local_vec(index);
f_a=f_local_vec(index+1);
u_z_max=-(chord(index+1)-chord(index-1))/2/(chord(index+1)-2*chord(index)+chord(index-1));
f_r=f_b+d_f/2+u_z_max*d_f;
w_r=2*pi*f_r;

%Obtain theta_r by linear interpolation
beta=atan2(imag(alpha_local_vec)-y_center,real(alpha_local_vec)-x_center);  %Returns the radial angles wrt the x-axis. Angles in the range [-pi~pi]
beta=unwrap(beta);  %corrects the phase angles (in radian) by adding multiples of 2pi when absolute jumps between consecutive elements are greater than or equal to the default jump tolerance

beta_b=beta(index);
beta_a=beta(index+1);
beta_r=beta_b+(f_r-f_b)*(beta_a-beta_b)/(f_a-f_b);

%Calculate eta_r (modal damping)
    %Method 1: average the two points before and after (for the solving in the final exam)
    w_b=2*pi*f_b;
    w_a=2*pi*f_a;
    theta_b=beta_b-beta_r;  %+ve
    theta_a=beta_r-beta_a;  %+ve
    eta_r=(w_a^2-w_b^2)/w_r^2/(tan(theta_a/2)+tan(theta_b/2));

    %Method 2: average points for which -pi/2 <= theta <= pi/2(for the solving in the final exam)
    %More accurate but longer method (Do not use this method in the final exam)
    theta_temp=beta-beta_r; %+ve direction is ccw, same as beta and opposite to the omega direction
    flag_b_vec=(theta_temp > 0) & (theta_temp <= pi/2);
    flag_a_vec=(theta_temp >= -pi/2) & (theta_temp < 0);

    theta_b_vec=theta_temp(flag_b_vec);  %+ve
    theta_a_vec=-theta_temp(flag_a_vec);  %+ve
    if length(theta_b_vec)<=1 || length(theta_a_vec)<=1,
        warning('Resolution of alpha_local_vec is coarse.'),
    end
    [theta_b_mat,theta_a_mat]=meshgrid(theta_b_vec,theta_a_vec);

    w_b_vec=2*pi*f_local_vec(flag_b_vec);
    w_a_vec=2*pi*f_local_vec(flag_a_vec);
    [w_b_mat,w_a_mat]=meshgrid(w_b_vec,w_a_vec);

    eta_r_mat=(w_a_mat.^2-w_b_mat.^2)/w_r^2./(tan(theta_a_mat/2)+tan(theta_b_mat/2));
    %eta_r=mean(mean(eta_r_mat));

% Mode shape calculation (A_r)
A_r=2*Radius*w_r^2*eta_r*exp(1i*(pi/2+beta_r));

% Residue calculation (B_r)
B_r=complex(x_center,y_center)-Radius*exp(1i*beta_r);

if ShowInternalDetails
    %For manual Display Only
    f_b
    f_a
    f_interp=f_local_vec(index-1:index+2).'
    f_mean=f_interp(1:end-1)+d_f/2
    z=chord(index-1:index+1)
    u_z_max
    f_r
    beta_b_rad=beta_b
    beta_b_deg=beta_b/pi*180
    beta_a_rad=beta_a
    beta_a_deg=beta_a/pi*180
    beta_r_rad=beta_r
    beta_r_deg=beta_r/pi*180
    theta_b_rad=theta_b
    theta_b_deg=theta_b/pi*180
    theta_a_rad=theta_a
    theta_a_deg=theta_a/pi*180
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Receptance_gen_local=A_r./(complex(w_r^2-(2*pi*f_local_vec).^2,eta_r*w_r^2))+B_r;
    CircleFitInternalDetails(f_local_vec,alpha_local_vec,Receptance_gen_local,circ_prop,f_r,beta_r,w_b_mat/2/pi,w_a_mat/2/pi,eta_r_mat,label_str);
end
