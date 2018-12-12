function [f_r,eta_r,A_r,B_r,circ_prop]=FRF_CircleFit(f_local_vec,Receptance_local_vec,ShowInternalDetails)

%Error checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the case of insufficient number of samples
N_pts = length(f_local_vec); 
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
[x_center,y_center,Radius] = err_fit_circle(real(Receptance_local_vec),imag(Receptance_local_vec));
circ_prop = struct('x_center',x_center,'y_center',y_center,'Radius',Radius);

%Identify the two points surrounding w_r
chord=sqrt((real(Receptance_local_vec(2:end))-real(Receptance_local_vec(1:end-1))).^2+(imag(Receptance_local_vec(2:end))-imag(Receptance_local_vec(1:end-1))).^2);
[chord_max,index]=max(chord);
f_before=f_local_vec(index);
f_after=f_local_vec(index+1);
u_z_max=-(chord(index+1)-chord(index-1))/2/(chord(index+1)-2*chord(index)+chord(index-1));
f_r=f_before+d_f/2+u_z_max*d_f;
w_r=2*pi*f_r;

%Obtain theta_r by linear interpolation
beta = atan2(imag(Receptance_local_vec)-y_center,real(Receptance_local_vec)-x_center);  %Returns the radial angles wrt the x-axis. Angles in the range [-pi~pi]
beta = unwrap(beta);  %corrects the phase angles (in radian) by adding multiples of 2pi when absolute jumps between consecutive elements are greater than or equal to the default jump tolerance

beta_before=beta(index);
beta_after=beta(index+1);
beta_r=beta_before+(f_r-f_before)*(beta_after-beta_before)/(f_after-f_before);

%Calculate eta_r (modal damping)
    %Fast method (for the solving in the final exam)
    w_before=2*pi*f_before;
    w_after=2*pi*f_after;
    theta_b=beta_before-beta_r;  %+ve
    theta_a=beta_r-beta_after;  %+ve
    eta_r=(w_after^2-w_before^2)/w_r^2/(tan(theta_a/2)+tan(theta_b/2));

    %More accurate but longer method (Do not use this method in the final exam)
    theta_temp=beta-beta_r; %+ve direction is ccw, same as beta and opposite to the omega direction
    flag_before= (theta_temp > 0) & (theta_temp <= pi/2);
    flag_after= (theta_temp >= -pi/2) & (theta_temp < 0);

    theta_before=theta_temp(flag_before);  %+ve
    theta_after=-theta_temp(flag_after);  %+ve
    [theta_before_mat,theta_after_mat]=meshgrid(theta_before,theta_after);

    w_before = 2*pi*f_local_vec(flag_before);
    w_after = 2*pi*f_local_vec(flag_after);
    [w_before_mat,w_after_mat]=meshgrid(w_before,w_after);

    eta_r_mat = (w_after_mat.^2-w_before_mat.^2)/w_r^2./(tan(theta_after_mat/2)+tan(theta_before_mat/2));

    %eta_r=mean(mean(eta_r_mat));

% Mode shape calculation (A_r)
A_r=2*Radius*w_r^2*eta_r*exp(1i*(pi/2+beta_r));

% Residue calculation (B_r)
B_r=complex(x_center,y_center)-Radius*exp(1i*beta_r);

if ShowInternalDetails
    %For manual Display Only
    f_before
    f_after
    f_interp=f_local_vec(index-1:index+2).'
    f_mean=f_interp(1:end-1)+d_f/2
    z=chord(index-1:index+1)
    u_z_max
    f_r
    beta_before_rad=beta_before
    beta_before_deg=beta_before/pi*180
    beta_after_rad=beta_after
    beta_after_deg=beta_after/pi*180
    beta_r_rad=beta_r
    beta_r_deg=beta_r/pi*180
    theta_b_rad=theta_b
    theta_b_deg=theta_b/pi*180
    theta_a_rad=theta_a
    theta_a_deg=theta_a/pi*180
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Receptance_gen_local=A_r./(complex(w_r^2-(2*pi*f_local_vec).^2,eta_r*w_r^2))+B_r;
    CircleFitInternalDetails(f_local_vec,Receptance_local_vec,Receptance_gen_local,circ_prop,beta_r,w_before_mat/2/pi,w_after_mat/2/pi,eta_r_mat);
end