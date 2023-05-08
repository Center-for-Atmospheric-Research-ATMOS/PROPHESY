# [336ed68f] CSV v0.10.9  from CSV v0.10.4 (if CSV version<v0.10.9, use string.(DataFrame()) instead of DataFrame())
# [a93c6f00] DataFrames v1.3.6

CSV.write(string(save_folder,"depth.csv"),DataFrame(r',:auto);header=true)
CSV.write(string(save_folder,"concentration_profile.csv"),DataFrame(ρA_1',:auto);header=true)
CSV.write(string(save_folder,"depth_lowres.csv"),DataFrame(r_lowres',:auto);header=true)

CSV.write(string(save_folder,"concentration_estimates.csv"),DataFrame(ρ_cp,:auto);header=true)
CSV.write(string(save_folder,"concentration_noise_variability_mean.csv"),DataFrame(mean_ρ_H',:auto);header=true)
CSV.write(string(save_folder,"concentration_noise_variability_cov.csv"),DataFrame(var_ρ_H,:auto);header=true)
CSV.write(string(save_folder,"concentration_distribution_mean_conditional_to_model.csv"),DataFrame(μρ_H',:auto);header=true)
CSV.write(string(save_folder,"concentration_distribution_cov_conditional_to_model.csv"),DataFrame(Γρ_H,:auto);header=true)

CSV.write(string(save_folder,"concentration_distribution_mean_one_model.csv"),DataFrame(μρ_HI',:auto);header=true)
for sample_flag in 1:Nsample
    CSV.write(string(save_folder,"concentration_distribution_mean_one_model_sample_",sample_flag,".csv"),DataFrame(Γρ_HI[:,:,sample_flag],:auto);header=true)
end
CSV.write(string(save_folder,"concentration_distribution_relative_energy_evolution_one_model.csv"),DataFrame(cumsum(-deltaU,dims=1)[1:1000:end,:],:auto);header=true)



CSV.write(string(save_folder,"concentration_model_variability_estimates.csv"),DataFrame(ρ_cp_HI,:auto);header=true)
CSV.write(string(save_folder,"concentration_model_variability_mean.csv"),DataFrame(mean_ρ_y',:auto);header=true)
CSV.write(string(save_folder,"concentration_model_variability_cov.csv"),DataFrame(var_ρ_y,:auto);header=true)
CSV.write(string(save_folder,"concentration_distribution_mean_conditional_to_noise.csv"),DataFrame(μρ_y',:auto);header=true)
CSV.write(string(save_folder,"concentration_distribution_cov_conditional_to_noise.csv"),DataFrame(Γρ_y,:auto);header=true)

CSV.write(string(save_folder,"concentration_distribution_mean_one_data.csv"),DataFrame(μρ_HI_sample',:auto);header=true)
for sample_flag in 1:Nsample
    CSV.write(string(save_folder,"concentration_distribution_mean_one_data_sample_",sample_flag,".csv"),DataFrame(Γρ_HI_sample[:,:,sample_flag],:auto);header=true)
end
CSV.write(string(save_folder,"concentration_distribution_relative_energy_evolution_one_data.csv"),DataFrame(cumsum(-deltaUh,dims=1)[1:1000:end,:],:auto);header=true)