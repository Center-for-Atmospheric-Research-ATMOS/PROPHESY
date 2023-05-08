function dataAndFit_xlsx2df(fileName::String;
    regData::Regex=r"Eph ?= ?",regFit::Regex=r"[Ff]itt?ing(_| )?results?",regEph::Regex=r"=[0-9]*eV",sortBy::Symbol=Symbol("Photon energy"),thenSortBy::Symbol=Symbol("Binding energy"))

    local xf_data = XLSX.readxlsx(fileName);
    local xf_data_sheet_names = XLSX.sheetnames(xf_data);
    local list_match = match.(regData,xf_data_sheet_names); 
    local list_match_fit = match.(regFit,xf_data_sheet_names);
    local xf_raw_data_sheet_names = xf_data_sheet_names[list_match.!=nothing];
    local xf_fit_data_sheet_names = xf_data_sheet_names[list_match_fit.!=nothing];

    Ndata = length(xf_raw_data_sheet_names);
    dictAllData = Dict();
    for xf_name in xf_raw_data_sheet_names
        # local x = XLSX.getdata(xf_data[xf_name])[2:end,:]
        # local col_sym = string.(XLSX.gettable(xf_data[xf_name])[2]) # valid for XLSX@v0.7.10 and under
        # local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,length(col_sym));
        # for j in eachindex(col_sym)
        #     if (typeof(x[1,j])<:AbstractString)
        #         dataPairs[j] = (col_sym[j] => parse.(Cdouble,x[:,j]))
        #     else
        #         dataPairs[j] = (col_sym[j] => convert(Array{Cdouble,1},x[:,j]))
        #     end
        # end
        local G = DataFrame(XLSX.gettable(xf_data[xf_name]))
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,ncol(G));
        for j in 1:ncol(G)
            if ((typeof(G[!,Symbol(names(G)[j])][1]))<:AbstractString)
                dataPairs[j] = (names(G)[j] => parse.(Cdouble,G[!,Symbol(names(G)[j])]))
            else
                dataPairs[j] = (names(G)[j] => convert(Array{Cdouble,1},G[!,Symbol(names(G)[j])]))
            end
        end
        df = DataFrame(dataPairs);
        dictAllData[Symbol(string("hν_",match.(regEph,xf_name).match[2:end-2]))] = df
    end


    x_fit = XLSX.getdata(xf_data[xf_fit_data_sheet_names[1]]);
    # remove missing columns
    x_fit = x_fit[:,broadcast(~,(ismissing.(x_fit[1, :])))]
    df_fit = DataFrame([col for col in eachcol(x_fit[2:end, :])], Symbol.(x_fit[1, :]))
    # remove missing rows
    # df_fit = dropmissing(df_fit, disallowmissing=true)
    filter!(x -> any(!ismissing, x), df_fit)
    df_fit = df_fit |> @orderby(_[sortBy]) |> @thenby(_[thenSortBy]) |> DataFrame

    # return
    dictAllData,df_fit,Ndata
end

function dataAndFit_xlsx2df_missing(fileName::String;
    regData::Regex=r"Eph ?= ?",regFit::Regex=r"[Ff]itt?ing(_| )?results?",regEph::Regex=r"=[0-9]*eV",sortBy::Symbol=Symbol("Photon energy"),thenSortBy::Symbol=Symbol("Binding energy"))

    local xf_data = XLSX.readxlsx(fileName);
    local xf_data_sheet_names = XLSX.sheetnames(xf_data);
    local list_match = match.(regData,xf_data_sheet_names); 
    local list_match_fit = match.(regFit,xf_data_sheet_names);
    local xf_raw_data_sheet_names = xf_data_sheet_names[list_match.!=nothing];
    local xf_fit_data_sheet_names = xf_data_sheet_names[list_match_fit.!=nothing];

    Ndata = length(xf_raw_data_sheet_names);
    dictAllData = Dict();
    for xf_name in xf_raw_data_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])[2:end,:]
        local col_sym = string.(XLSX.gettable(xf_data[xf_name])[2])
        local dataPairs = Array{Pair{String,Vector{Union{Missing, Float64}}}}(undef,length(col_sym))
        for j in 1:length(col_sym)
            local y = Array{Union{Missing, Float64}, 1}(undef,length(x[:,j]));
            for i in 1:length(x[:,j])
                if (typeof(x[i,j])!=Missing)
                    if (typeof(x[i,j])<:AbstractString)
                        y[i] = parse(Cdouble,x[i,j])
                    else
                        y[i] = convert(Cdouble,x[i,j])
                    end
                else
                    y[i] = missing
                end
            end
            dataPairs[j] = (col_sym[j] => y)
        end
        df = DataFrame(dataPairs);
        dictAllData[Symbol(string("hν_",match.(regEph,xf_name).match[2:end-2]))] = df
    end


    x_fit = XLSX.getdata(xf_data[xf_fit_data_sheet_names[1]]);
    # remove missing columns
    x_fit = x_fit[:,broadcast(~,(ismissing.(x_fit[1, :])))]
    df_fit = DataFrame([col for col in eachcol(x_fit[2:end, :])], Symbol.(x_fit[1, :]))
    # remove missing rows
    # df_fit = dropmissing(df_fit, disallowmissing=true)
    filter!(x -> any(!ismissing, x), df_fit)
    df_fit = df_fit |> @orderby(_[sortBy]) |> @thenby(_[thenSortBy]) |> DataFrame

    # return
    dictAllData,df_fit,Ndata
end

function curveFromFit(dfPeak::DataFrame,Be::Array{Cdouble,1},shiftPeaks::Bool;
    bind_sym::Symbol=Symbol("Binding energy"), shift_sym::Symbol=Symbol("Peak shift"), gauss_sym::Symbol=Symbol("FWHM(G)"), area_sym::Symbol=Symbol("Area"), loren_sym::Symbol=Symbol("FWHM(L)"))
    local Npeaks   = length(dfPeak[!,bind_sym]);
    local BePeak   = zeros(Cdouble,Npeaks);
    local σePeak   = zeros(Cdouble,Npeaks);
    local σePeakL  = zeros(Cdouble,Npeaks);
    local AePeak   = zeros(Cdouble,Npeaks);
    local σ_peak   = zeros(Cdouble,Npeaks,length(Be));

    for k in 1:Npeaks
        if (typeof(dfPeak[!,bind_sym][k])<:AbstractString)
            BePeak[k] = parse(Cdouble,dfPeak[!,bind_sym][k])
        else
            BePeak[k] = dfPeak[!,bind_sym][k]
        end
        if shiftPeaks
            if (typeof(dfPeak[!,shift_sym][k])<:AbstractString)
                BePeak[k] = BePeak[k] + parse(Cdouble,dfPeak[!,shift_sym][k])
            else
                BePeak[k] = BePeak[k] + dfPeak[!,shift_sym][k]
            end
        end
        if (typeof(dfPeak[!,gauss_sym][k])<:AbstractString)
            σePeak[k] = 0.5*1.0e-3parse(Cdouble,dfPeak[!,gauss_sym][k])
        else
            σePeak[k] = 0.5*1.0e-3dfPeak[!,gauss_sym][k];
        end
        if (typeof(dfPeak[!,area_sym][k])<:AbstractString)
            AePeak[k] = parse(Cdouble,dfPeak[!,area_sym][k])
        else
            AePeak[k] = dfPeak[!,area_sym][k]
        end
        if (typeof(dfPeak[!,loren_sym][k])<:AbstractString)
            σePeakL[k] = 0.5*1.0e-3parse(Cdouble,dfPeak[!,loren_sym][k])
        else
            σePeakL[k] = 0.5*1.0e-3dfPeak[!,loren_sym][k];
        end
        σ_peak[k,:] = (1.0/sqrt(2π*σePeak[k]^2))*exp.(-0.5*((Be.-BePeak[k])/σePeak[k]).^2)
    end
    # return
    Npeaks, BePeak, σePeak, σePeakL, AePeak, σ_peak
end

















function dataAndMeta_xlsx2df(fileName::String,regData::Regex,regMeta::Regex)
    local xf_data = XLSX.readxlsx(fileName);
    local xf_sheet_names = XLSX.sheetnames(xf_data);
    local list_match_data = match.(regData,xf_sheet_names);
    local list_match_meta = match.(regMeta,xf_sheet_names);
    local xf_data_sheet_names = xf_sheet_names[list_match_data.!=nothing];
    local xf_meta_sheet_names = xf_sheet_names[list_match_meta.!=nothing];
    # the data
    dictAllData = Dict();
    symbolDict = Dict();
    for xf_name in xf_data_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,size(x,2));
        for j in 1:size(x,2)
            if (typeof(x[2:end,j][1])<:AbstractString)
                local data_x = coalesce.(x[2:end,j],"0.0") # [.!ismissing.(x[2:end,j])]; # eltype(x)
                dataPairs[j] = (x[1,j] => parse.(Cdouble,replace.(data_x, r"[,;]" => "")))
            else
                local data_x = coalesce.(x[2:end,j],0.0)
                dataPairs[j] = (x[1,j] => convert(Array{Cdouble,1},data_x))
            end
        end
        df = DataFrame(dataPairs);
        dictAllData[Symbol(xf_name)] = df 
        symbolDict[Symbol(xf_name)] = Symbol(string("λe_",df[!,:λ][1]))
    end
    # meta data
    local dictAllMeta = Dict();
    for xf_name in xf_meta_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,size(x,2));
        for j in 1:size(x,2)
            dataPairs[j] = (x[1,j] => convert(Array{Cdouble,1},x[2:end,j]))
        end
        df = DataFrame(dataPairs);
        dictAllMeta[Symbol(xf_name)] = df 
    end

    local dfMeta = DataFrame();
    for (_,df) in dictAllMeta
        append!(dfMeta,df);
    end
    df = dfMeta |> @orderby(_[:hν]) |> DataFrame

    # return
    dictAllData,df,symbolDict
end

function model_xlsx2df(fileName::String)
    local xf_model = XLSX.readxlsx(fileName);
    local xf_sheet_names = XLSX.sheetnames(xf_model);
    # the data
    dictAllGeom = Dict();
    symbolDict = Dict();
    for xf_name in xf_sheet_names
        local x = XLSX.getdata(xf_model[xf_name])
        local dataPairs = Array{Pair{String,Union{Vector{Cdouble},Vector{String}}}}(undef,size(x,2));
        for j in 1:size(x,2)
            if (x[1,j]=="model")
                dataPairs[j] = (x[1,j] => string.(x[2:end,j]))
            else
                dataPairs[j] = (x[1,j] => convert(Array{Cdouble,1},x[2:end,j]))
            end
        end
        df = DataFrame(dataPairs);
        dictAllGeom[Symbol(xf_name)] = df 
        symbolDict[Symbol(xf_name)] = Symbol(string("hν_",round(Int64,df[!,:hν][1])))
    end

    # return
    dictAllGeom,symbolDict
end


function IGORcolumn_xlsx2df(fileName::String,regSheet::Regex,regData::Regex,regMeta::Regex)
    local xf_data = XLSX.readxlsx(fileName);
    local xf_sheet_names = XLSX.sheetnames(xf_data);
    local list_match_data = match.(regSheet,xf_sheet_names);
    local xf_data_sheet_names = xf_sheet_names[list_match_data.!=nothing];

    # the data
    dictAllData = Dict();
    dictAllMeta = Dict();
    for xf_name in xf_data_sheet_names
        local x = XLSX.getdata(xf_data[xf_name])

        ## data
        local list_data_column = match.(regData,coalesce.(x[1,:],""));
        local x_data = x[:,list_data_column.!=nothing];
        local dataPairs = Array{Pair{String,Vector{Cdouble}}}(undef,size(x_data,2));
        for j in 1:size(x_data,2)
            if (typeof(x_data[2:end,j][1])<:AbstractString)
                local data_x = coalesce.(x_data[2:end,j],"0.0") # [.!ismissing.(x[2:end,j])]; # eltype(x)
                dataPairs[j] = (x_data[1,j] => parse.(Cdouble,replace.(data_x, r"[,;]" => "")))
            else
                local data_x = coalesce.(x_data[2:end,j],0.0)
                dataPairs[j] = (x_data[1,j] => convert(Array{Cdouble,1},data_x))
            end
        end
        df = DataFrame(dataPairs);
        dictAllData[Symbol(xf_name[1:end-5])] = df 

        ## meta data
        local list_meta_column = match.(regMeta,coalesce.(x[1,:],""));
        local x_meta = x[:,list_meta_column.!=nothing];
        local metaPairs = Array{Pair{String,Vector{Cdouble}}}(undef,size(x_meta,2));
        for j in 1:size(x_meta,2)
            meta_x = collect(skipmissing(x_meta[2:end,j]))
            if (typeof(meta_x[1])<:AbstractString)
                # local data_x = coalesce.(x_meta[2:end,j],"0.0") # [.!ismissing.(x[2:end,j])]; # eltype(x)
                metaPairs[j] = (x_meta[1,j] => parse.(Cdouble,replace.(meta_x, r"[,;]" => "")))
            else
                # local data_x = coalesce.(x_meta[2:end,j],0.0)
                metaPairs[j] = (x_meta[1,j] => convert(Array{Cdouble,1},meta_x))
            end
        end
        dictAllMeta[Symbol(xf_name[1:end-5])] = DataFrame(metaPairs);
    end

    # return
    dictAllData,dictAllMeta
end