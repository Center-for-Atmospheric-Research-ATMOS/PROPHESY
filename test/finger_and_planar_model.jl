
H_z_far,H_zxy_far,_,_,_   = plane_gain_H(x_far,Y,z_surf,x0_far,y0_far,z0_far,λe);
H_z_near,H_zxy_near,_,_,_ = plane_gain_H(x_near,Y,z_surf,x0_near,y0_near,z0_near,λe);
H_z_finger_far,_          = finger_gain_H(0.0,0.0,z_surf,0.0,0.0,z0_far,λe);
H_z_finger_near,_         = finger_gain_H(0.0,0.0,z_surf,0.0,0.0,z0_near,λe);
# H_z_z_0,_                 = finger_gain_H(0.0,0.0,z_surf,0.0,0.0,z0_far,0.58λe);

relative_eal_far  = log.(H_z_finger_far[2:end-1]/H_z_finger_far[end-1])./log.(H_z_far[2:end-1]/H_z_far[end-1]);
relative_eal_near = log.(H_z_finger_near[2:end-1]/H_z_finger_near[end-1])./log.(H_z_near[2:end-1]/H_z_near[end-1]);
r_eal_min_far,r_eal_max_far   = extrema(relative_eal_far[1:end-1])
r_eal_min_near,r_eal_max_near = extrema(relative_eal_near[1:end-1])

H_z_far_min,_          = finger_gain_H(x0_far,y0_far,z_surf,x0_far,y0_far,z0_far,r_eal_min_far*λe);
H_z_far_max,_          = finger_gain_H(x0_far,y0_far,z_surf,x0_far,y0_far,z0_far,r_eal_max_far*λe);
H_z_near_min,_         = finger_gain_H(x0_near,y0_near,z_surf,x0_near,y0_near,z0_near,r_eal_min_near*λe);
H_z_near_max,_         = finger_gain_H(x0_near,y0_near,z_surf,x0_near,y0_near,z0_near,r_eal_max_near*λe);


figure();
plot(z_surf[2:end-1],relative_eal_far, label="far")
plot(z_surf[2:end-1],relative_eal_near, label="near")
legend()
xlabel("depth [\$\\mu m\$]")
ylabel("log ratio: \$\\frac{\\log(\\frac{I^{\\mathrm{planar}}}{I^{\\mathrm{planar}}_0} )}{\\log(\\frac{I^{\\mathrm{cylindrical}}}{I^{\\mathrm{cylindrical}}_0})}\$")

## plot
figure();
plot(z_surf,H_z_far/maximum(H_z_far), label="planar: far \$\\lambda_e\$", color="tab:blue");
plot(z_surf,H_z_near/maximum(H_z_near), label="planar: near \$\\lambda_e\$", color="tab:green");
plot(z_surf,H_z_finger_far/maximum(H_z_finger_far), label="finger approximation: far \$\\lambda_e\$", color="tab:orange");
plot(z_surf,H_z_finger_near/maximum(H_z_finger_near), label="finger approximation: near \$\\lambda_e\$", color="tab:olive");
s = @sprintf "finger approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
fill_between(z_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:red",label=s)
s = @sprintf "finger approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
fill_between(z_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:pink",label=s)
legend()
xlabel("depth [\$\\mu m\$]")
ylabel("normalized gain")
