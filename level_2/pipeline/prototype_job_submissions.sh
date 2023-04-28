SEED=`cat gt_seed`

julia network_estimation.jl *.csv *.tre 0 3 `expr $SEED + 1` > snaq_outgroup_h0.outerr 2>&1 

julia network_estimation.jl *.csv snaq_output_h0.out 1 3 20 `expr $SEED + 3` > snaq_outgroup_h1.outerr 2>&1

julia network_estimation.jl *.csv snaq_output_h1.out 2 3 20 `expr $SEED + 5` > snaq_outgroup_h2.outerr 2>&1

julia network_estimation.jl *.csv snaq_output_h2.out 3 3 20 `expr $SEED + 7` > snaq_outgroup_h3.outerr 2>&1 
