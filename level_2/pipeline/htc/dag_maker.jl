parent_dir = "../output/all/jobs"

missing_job_nos = Vector{Int}()

for subdir_no in 0:162000
    if !isfile("$parent_dir/job_$(subdir_no)/emp_est_out$(subdir_no).tar.gz")
        push!(missing_job_nos, subdir_no)
    end
end

datasets = missing_job_nos[1:10000]
#datasets = [0,1]
h_steps = [0, 1, 2, 3, 4, 5]

# Open the file for writing
open("./htc/master_workflow.dag", "w") do f
    write(f, "CONFIG dagman.config\n")
    write(f, "MAXJOBS PREP_PHASE 200\n")
    write(f,"RETRY ALL_NODES 2\n\n")
    for data in datasets
        
        # --- Define the CF Prep Job ---
        prep_name = "$(data)_prep"
        write(f, "JOB  $prep_name  gen_cfs.sub\n")
        write(f, "VARS $prep_name  job_id=\"$data\"\n\n")
        write(f, "CATEGORY $prep_name PREP_PHASE\n\n")
        
        # Set the prep job as the parent for the h=0 step!
        previous_job = prep_name 
        # -----------------------------------
        
        for h in h_steps
            job_name = "$(data)_h$(h)"
            stats_name = "$(data)_h$(h)_stats"

            job_priority = h+1
            
            # Route the start tree
            if h == 0
                start_tree = "qmc_tree.newick"
            else
                start_tree = "h_$(h-1).out"
            end
            
            # Write JOB definition
            write(f, "JOB  $job_name  snaq_runs.sub\n")
            write(f, "VARS $job_name  h_val=\"$h\" job_id=\"$data\" start_tree=\"$start_tree\"\n")
            write(f, "PRIORITY $job_name $job_priority\n")
            
            # POST script arguments
            write(f, "SCRIPT POST $job_name  /home/justison/scripts/shells/find_best.sh $h $data\n")
            write(f, "PARENT $previous_job CHILD $job_name\n")
            

            write(f, "\n")
            previous_job = job_name
        end

        # --- The Final Stats Job ---
        stats_name = "$(data)_final_stats"
        write(f, "JOB  $stats_name  snaq_analysis.sub\n")
        write(f, "VARS $stats_name  job_id=\"$data\"\n")
        write(f, "PRIORITY $stats_name 10\n")
        
        # PRE script to package the networks
        write(f, "SCRIPT PRE $stats_name /home/justison/scripts/shells/compress_nets.sh $data\n")

        # Link it to the very last h-step!
        write(f, "PARENT $previous_job CHILD $stats_name\n\n")
        write(f, "SCRIPT POST $stats_name  /home/justison/scripts/shells/post_bookkeeping.sh $data\n")

    end
end

println("DAG file created: master_workflow.dag for $(length(datasets)) missing jobs.")