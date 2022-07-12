function preparehdf(hdfname_i::String, hdfname_ii::String, n_cells::Int, n_el::Int, n_gp::Int)
    εⁱ_0 = zeros(6,8*n_cells) # abh. Größe RVE nElmt*nGP - Anfangsbelegung = 0
    work_dir = pwd()
    cd("Results")
    h5open(hdfname_i, "w") do file
        for i in 1:(n_el*n_gp)
            dataname = "epsi_" * string(i)
            write(file, dataname, εⁱ_0)
        end         
    end
    h5open(hdfname_ii, "w") do file
        for i in 1:(n_el*n_gp)
            dataname = "epsi_" * string(i)
            write(file, dataname, εⁱ_0)
        end         
    end
    cd(work_dir) 
end

function overwritehdf(hdfname::String, dataname::String, data::Array{Float64, 2})
    work_dir = pwd()
    cd("Results")
    h5open(hdfname, "r+") do file 
        dataset = file[dataname]
        dataset[:,:] = data[:,:]
    end
    cd(work_dir)
    return
end

function readhdf(hdfname::String, dataname::String)
    work_dir = pwd()
    cd("Results")
    data = h5open(hdfname, "r") do file
        read(file, dataname)
    end
    cd(work_dir)
    return data
end

function updatehdf(hdf_i::String, hdf_ii::String) # delete _i file, copy ii to i 
    work_dir = pwd()
    cd("Results")
    if(isfile(hdf_i))
        rm(hdf_i)
    end
    cp(hdf_ii, hdf_i)
    cd(work_dir)
    return 
end