# write get-pop-tree() for julia
function get_pop_tree(csv_file, qmc_file, tre_file, columns)
    # run qmcformat to create the proper file
    qmcformat(csv_file, qmc_file, columns, "qmc_file")
    
    # pick the right binary
    if Sys.isapple()
        qmc_bin = "find-cut-Mac"
    elseif Sys.islinux()
        qmc_bin = "find-cut-Linux-64"
    else
        error("QuartetMaxCut only supports either Linux or Mac operating systems")
    end
    # run the program on the input qmc file and return to the output tree file
    run(`$qmc_bin qrtt="$qmc_file" otre="$tre_file"`)
    println("Finished running QMC")
    # reformat QMC's output in order to back-translate the original leaf labels
    dict_taxa = qmcformat(csv_file, qmc_file, columns, "dict")
    # reverse the dict so that the indices are now keys and labels are values
    #dict_taxa = Dict(string.(values(dict_taxa)) .=> keys(dict_taxa))
    # dicts are not necessarily the best option here but they are already working up the line, ket the values for yeis and idndices
    dict_keys = collect(keys(dict_taxa))
    dict_values = collect(values(dict_taxa))
    # read the tre file just generated
    tree = readline(tre_file)
    println(tree)
    # iterate over keys and do replacements
    for i in 1:length(dict_keys)
        println(i)
        tree = replace(tree, string(dict_values[i]) => string(dict_keys[i]))
    end
    println(tree)
end

# test get_pop_tree
get_pop_tree("test.csv", "test.qmc", "test.tre", [1, 2, 3, 4, 5, 8, 11])
get_pop_tree("1_seqgen.CFs.csv", "1_seqgen.CFs.qmc", "1_seqgen.CFs.tre", [1, 2, 3, 4, 5, 8, 11])
