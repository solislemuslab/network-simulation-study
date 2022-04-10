# SNaQ running time study

Goal: Study the variables that play a role in the running time of SNaQ.

## Data
- We will use the simulated data for snaq and ransanec projects (`knownGT`)
- We have 5 networks of increasing number of tips: n6, n10, n15, n25, n50
- For each network (n6,n10,n15), we have 10 replicates obtained as follow:
    - 2 replicates per number of gene trees: 30,100,300,1000,3000 (`1_astral.in` and `2_astral.in` renamed as `1_astral_30.in` for 30 gene trees e.g.)
    - we do not expect to see changes in running time due to number of gene trees
- For each network (n25, n50), we have 10 replicates obtained as follow:
    - 2 replicates per number of gene trees: 50, 100, 500, 1000, 5000 e.g. `n25h5-gt50-1.tre`, `n25h5-gt50-2.tre`
    - Note that files with 5000 gene trees had to be saved into 5 files of 1000 gene trees (see more details below).


## Running SNaQ

- Maybe we want to run it on chtc?
    - GAB: I'm not sure, maybe the runtime will be affected by processor architecture, and the way CHTC is built we cannot be sure of whether individual runs went to different processors and thus are showing variation because of this. My feeling is that we need to run all of them on the same machine, or a couple times in different machines so that we can be sure that under the same architecture we are really characterising the runtime as a function of our parameters of interest.

1. Create CF tables (save running time too and memory also)
2. Run snaq in parallel for h=0 and nruns=1

### Temp trials before running on WID server

We want to use the julia scripts:
- `cfs_and_starting_tree.jl`
- `network_estimation.jl`
in `snaq_outgroup` so that we do not have repeated identical functions.

We had to make sure that `get-pop-tree.pl` was an executable and it was in the path, and that QMC was also in the path.

We go into the `running-time/data/n6` folder and do one small trial:
```shell
julia ../../../snaq_outgroup/cfs_and_starting_tree.jl 1_astral_30.in
```

This ran and created two files:
- `CFs_1_astral_30.in.csv`
- `CFs_1_astral_30.in.csv.QMC.tre`

Now, we need to know how to check time and memory for this command.

I had to install the GNU time OSX command with `brew install gnu-time` and it is used as `gtime`:
```shell
$ gtime -f "mem=%K RSS=%M elapsed=%E cpu.sys=%S user=%U" julia ../../../snaq_outgroup/cfs_and_starting_tree.jl 1_astral_30.in

Processing 1_astral_30.in...
Reading in trees, looking at 15 quartets in each...
0+------------------------------+100%
  ******************************
30.0 gene trees per 4-taxon set

Script was called as follows:
perl get-pop-tree.pl CFs_1_astral_30.in.csv

Parsing major resolution of each 4-taxon set... done.
Running Quartet Max Cut...

Quartet MaxCut version 2.10 by Sagi Snir, University of Haifa

quartet file is CFs_1_astral_30.in.csv.QMC.txt, 

Number of quartets is 15, max element 6

Number of quartets read: 15, max ele 6

Started working at  Sun Apr 10 08:37:59 2022

Ended working at  Sun Apr 10 08:37:59 2022

Quartet Max Cut complete, tree located in 'CFs_1_astral_30.in.csv.QMC.tre'.

mem=0 RSS=623016 elapsed=0:13.04 cpu.sys=0.56 user=13.09
```

We will search exactly what we want in the [time manual](https://man7.org/linux/man-pages/man1/time.1.html):

- %M     Maximum resident set size of the process during its lifetime in Kbytes.
- %K     Average total (data+stack+text) memory use of the process, in Kbytes.
- %E     Elapsed real time (in [hours:]minutes:seconds)
- %e     Elapsed real time (in seconds)
- %S     Total number of CPU-seconds that the process spent in kernel mode.

```shell
$ gtime -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S" julia ../../../snaq_outgroup/cfs_and_starting_tree.jl 1_astral_30.in
```

We want to make sure to send things to files, so:

```shell
$ gtime -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S" julia ../../../snaq_outgroup/cfs_and_starting_tree.jl 1_astral_30.in 2> 1_astral_30.time 1> 1_astral_30.log
```


**To do** We need to decide which memory options to keep (discussing with Gustavo and Carlos).

### Concatenating gene trees for n25 and n50

We used the [ransanec](https://github.com/crsl4/ransanec) simulated gene trees for the cases of `n25` and `n50`. 
There was a bug in the simulating script simulating-gene-trees.jl where the same set of gene trees was written over multiple output files. For some reason, I was unable to produce a file with over 1000 gene trees, so for 5000 gene trees, I had to produce 5 files. But because of the bug, all 5 files were identical.

Thus, we deleted all the gene trees files with 5000 trees (until we can fix those simulations). Another thing to note is that we do not have gene tree file for 1000 gene trees for `n50`, so we will use one of the files from 5000 (that has only 1000 gene trees):
- `n50h10-gt5000-1-1.tre` -> `n50h10-gt1000-1.tre`
- `n50h10-gt5000-2-1.tre` -> `n50h10-gt1000-2.tre`

**To do** Waiting for simulations for 5000 gene trees for n25 and n50, and concatenate them.

## Running time and memory for creation of CF tables

We will run these in my Mac. We need information about the machine (see [here](https://vitux.com/get-linux-system-and-hardware-details-on-the-command-line/)):

```shell
$ uname -a
Darwin C02Z60BVM0XV.local 20.6.0 Darwin Kernel Version 20.6.0: Mon Aug 30 06:12:21 PDT 2021; root:xnu-7195.141.6~3/RELEASE_X86_64 x86_64

$ sysctl -a | grep cpu | grep hw 
hw.ncpu: 36
hw.activecpu: 36
hw.cpu64bit_capable: 1
hw.cpufamily: 939270559
hw.cpufrequency: 2300000000
hw.cpufrequency_max: 2300000000
hw.cpufrequency_min: 2300000000
hw.cpusubfamily: 0
hw.cpusubtype: 8
hw.cputype: 7
hw.logicalcpu: 36
hw.logicalcpu_max: 36
hw.physicalcpu: 18
hw.physicalcpu_max: 18
hw.cputhreadtype: 1
```

Also, from the `About my Mac`:
```
Hardware Overview:

  Model Name:	iMac Pro
  Model Identifier:	iMacPro1,1
  Processor Name:	18-Core Intel Xeon W
  Processor Speed:	2.3 GHz
  Number of Processors:	1
  Total Number of Cores:	18
  L2 Cache (per Core):	1 MB
  L3 Cache:	24.8 MB
  Hyper-Threading Technology:	Enabled
  Memory:	256 GB
  System Firmware Version:	1554.140.20.0.0 (iBridge: 18.16.14759.0.1,0)
  Serial Number (system):	C02Z60BVM0XV
  Hardware UUID:	56C27A24-0BBC-587A-94B5-3A6058911698
  Provisioning UDID:	56C27A24-0BBC-587A-94B5-3A6058911698
  Activation Lock Status:	Enabled
```

Let's start with n6, n10, n15 that have similar structure.
We need to do three loops:
- n6, n10, n15
- 30, 100, 300, 1000, 3000
- 1, 2 (replicate)

```shell
for i in 6 10 15
do
echo "n$i"
done
```

We want to be inside the `running-time/data` folder. Testing:
```shell
for i in 6 10 15
do
for j in 30 100 300 1000 3000
do
for k in 1 2
do
echo "$i, $j, $k"
echo "julia ../../snaq_outgroup/cfs_and_starting_tree.jl n$i/${k}_astral_${j}.in 2> n$i/${k}_astral_${j}.time 1> n$i/${k}_astral_${j}.log"
done
done
done
```

Real runs (started 4/10 9:36am):
```shell
for i in 6 10 15
do
for j in 30 100 300 1000 3000
do
for k in 1 2
do
gtime -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S" julia ../../snaq_outgroup/cfs_and_starting_tree.jl n$i/${k}_astral_${j}.in 2> n$i/${k}_astral_${j}.time 1> n$i/${k}_astral_${j}.log
done
done
done
```

We get an error in all the cases:
```
LoadError: SystemError: opening file "CFs_n6/1_astral_300.in.csv": No such file or directory
```

It seems that we need to be inside each of the folders. So, we had to delete all log and time output files.

Real runs (started 4/10 9:45am):
```shell
for i in 6 10 15
do
cd n$i
for j in 30 100 300 1000 3000
do
for k in 1 2
do
echo "$i, $j, $k"
gtime -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S" julia ../../../snaq_outgroup/cfs_and_starting_tree.jl ${k}_astral_${j}.in 2> ${k}_astral_${j}.time 1> ${k}_astral_${j}.log
done
done
cd ..
done
```
(finished at 9:51am)

For the cases of n25 and n50, we first have to rename the files so that it is easier to call them:
- `n25h5-gt50-1` -> `n25-gt50-1`

Real runs (started 4/10 9:55am):
```shell
for i in 25 50
do
cd n$i
for j in 50 100 500 1000
do
for k in 1 2
do
echo "$i, $j, $k"
gtime -f "mem=%K RSS=%M elapsed=%e cpu.sys=%S" julia ../../../snaq_outgroup/cfs_and_starting_tree.jl n${i}-gt${j}-${k}.tre 2> n${i}-gt${j}-${k}.time 1> n${i}-gt${j}-${k}.log
done
done
cd ..
done
```
(finished 10:03am)

**To do** 
- test the julia script (save time and memory also) on each of the cases
- move everything to the server to run in parallel there