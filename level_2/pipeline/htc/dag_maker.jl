parent_dir = "../output/all/jobs"


function get_original_job_id(setting_no, phy_no, rep_no)
    s_idx = setting_no - 1
    p_idx = phy_no - 1
    r_idx = rep_no - 1
    
    return r_idx + (p_idx * 30) + (s_idx * 150 * 30)
end

ordered_job_ids = Int[]
sizehint!(ordered_job_ids, 18 * 150 * 30)

for r in 1:30          # Replicates change slowest
    for p in 1:150     # Then phylogenies
        for s in 1:18  # Settings change fastest
            push!(ordered_job_ids, get_original_job_id(s, p, r))
        end
    end
end

missing_job_nos = Vector{Int}()
for subdir_no in ordered_job_ids
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

### create jobs tar to send to hpc
# Assuming ordered_job_ids contains your integers (e.g., [0, 4500, 9000, ...])
# Construct the relative paths: "jobs/job_0", "jobs/job_4500", etc.
# 1. Get the existing paths
file_paths = ["../jobs/job_$id" for id in datasets]
existing_paths = filter(isdir, file_paths)

# 2. Format the manifest paths to look like "jobs/job_XXXXX"
# This tells tar to build a "jobs" folder inside the archive
job_folders = ["jobs/$(basename(path))" for path in existing_paths]

manifest_file = "tar_manifest.txt"
open(manifest_file, "w") do io
    for folder in job_folders
        println(io, folder)
    end
end

# 3. Change directory to ".." (the parent of jobs) before archiving
# This keeps the paths relative to the parent and perfectly clean!
run(`tar -czf subset_jobs.tar.gz -C .. -T $manifest_file`)

# 4. Clean up the temporary manifest file
rm(manifest_file)