# finger_and_sphere_model.jl

H_r_sphere_far,H_rφθ_sphere_far,_,_,_   = sphere_gain_H(r_surf,φ_sphere,θ_sphere,x0_far,0.0,z0_far,μ0,λe);
H_r_sphere_near,H_rφθ_sphere_near,_,_,_ = sphere_gain_H(r_surf,φ_sphere,θ_sphere,x0_near,0.0,z0_near,μ0,λe);