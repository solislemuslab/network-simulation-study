#!/usr/bin/bash

# # custom network, not triggering ultrametricity warning, labeling all nodes as suggested by Joe
#  "(((h1#0.4:0.4,D:0.4)n3:0.3,((A:0.2,B:0.2)n2:0.3,(C:0.2)h1#0.6:0.3)n1:0.2)n0:0.2,E:0.9)r;"
# # custom network, converted from standard extended newick to hl via PhyloNetworks.hybridlambdaformat
#  "(((H1#0.5:4.0,D:4.0)I1:3.0,((A:2.0,B:2.0)I2:3.0,(C:2.0)H1#0.5:3.0)I3:2.0)I4:2.0,E:9.0)I5;"
# # custom network, not triggering ultrametricity warning, labeling all nodes as suggested by Joe
#  "(((h1#0.5:4.0,D:4.0)n3:3.0,((A:2.0,B:2.0)n2:3.0,(C:2.0)h1#0.5:3.0)n1:2.0)n0:2.0,E:9.0)r;"
# # custom network, triggering ultrametricity warning, labeling all nodes as suggested by Joe
#  "(((h1#0.5:2.0,D:4.0)n3:3.0,((A:2.0,B:2.0)n2:3.0,(C:2.0)h1#0.5:3.0)n1:2.0)n0:2.0,E:9.0)r;"
# # custom network, not triggering ultrametricity warning, with or without root string
#  "(((h1#0.5:4.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0);"
#  "(((h1#0.5:4.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0)r;"
# # custom network, triggering ultrametricity warning, with or without root string
#  "(((h1#0.5:2.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0);"
#  "(((h1#0.5:2.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0)r;"
# # 3_tax_multi_tax_test_hybrid1, with or without root string
#  "(((B:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6)r;"
#  "(((B:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6);"
# # 6_tax_multi_tax_test_hybrid1_topo1, with or without root string
#  "((((B:1,(D:0.6,(E:0.6,F:0.6)s5:0.6)s4:0.6)s3:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6)r;"
#  "((((B:1,(D:0.6,(E:0.6,F:0.6)s5:0.6)s4:0.6)s3:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6);"

VARS=("(((h1#0.4:0.4,D:0.4)n3:0.3,((A:0.2,B:0.2)n2:0.3,(C:0.2)h1#0.6:0.3)n1:0.2)n0:0.2,E:0.9)r;" "(((H1#0.5:4.0,D:4.0)I1:3.0,((A:2.0,B:2.0)I2:3.0,(C:2.0)H1#0.5:3.0)I3:2.0)I4:2.0,E:9.0)I5;" "(((h1#0.5:4.0,D:4.0)n3:3.0,((A:2.0,B:2.0)n2:3.0,(C:2.0)h1#0.5:3.0)n1:2.0)n0:2.0,E:9.0)r;" "(((h1#0.5:2.0,D:4.0)n3:3.0,((A:2.0,B:2.0)n2:3.0,(C:2.0)h1#0.5:3.0)n1:2.0)n0:2.0,E:9.0)r;" "(((h1#0.5:4.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0);" "(((h1#0.5:4.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0)r;" "(((h1#0.5:2.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0);" "(((h1#0.5:2.0,D:4.0):3.0,((A:2.0,B:2.0):3.0,(C:2.0)h1#0.5:3.0):2.0):2.0,E:9.0)r;" "(((B:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6)r;" "(((B:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6);" "((((B:1,(D:0.6,(E:0.6,F:0.6)s5:0.6)s4:0.6)s3:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6)r;" "((((B:1,(D:0.6,(E:0.6,F:0.6)s5:0.6)s4:0.6)s3:.6)h1#.5:.6,A:1.2)s1:.6,(h1#.5:.6,C:1.2)s2:.6);")

for i in ${VARS[@]}; do
    echo $i
done

for i in ${VARS[@]}; do
    printf "\nProcessing $i...\n"
    hybrid-Lambdav0.6.2-pathdiff -spcu $i -num 3 -seed 2 -o out0.6.2
    grep "#" out0.6.2_coal_unit
#    cat out0.6.2_coal_unit
    hybrid-Lambdav0.6.3 -spcu $i -num 3 -seed 2 -o out0.6.3
    grep "#" out0.6.3_coal_unit
#    cat out0.6.3_coal_unit
done

rm out0*
