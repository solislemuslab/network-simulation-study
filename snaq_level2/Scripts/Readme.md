# Pipeline for assessing the behavior of snaq on level-2 networks

The tasks accomplished by each of the scripts is described below.. Please note that execution should be carried out as indicated.

## `Species_Gene_Networks.sh`

**What the script does.**
- Create all the path structures to store all the results in an orderly way.

- Simulate species Networks with SimSpeciesNetworks.r code.

- Select only the species Networks with no single reticulation (In the future we are going to consider all the species Networks).

- Save the spesie Networks in the temporary path "Output/RNetworks".

- Read the simulated species Networks with julia, save it as julia format in the temporary path "Output/RJuliaNetworks" and in the Hybrid lambda format in the temporary path "Output/JuliaForHybLam".

- Copy and paste the species Networks in Hybrid lambda format to "src" of hybrid lambda path, run Hybrid lambda and move all the results into the temporary path "Output/TemporaryHyLambdaStorage".

- Evaluate the Hybrid lambda results to know if the specie networks are considered for Hybrid lambda as Ultrametric or not and save the previous results in the path "Output/UltraMetric" and "Output/NonUltraMetric" as appropriate using the script UltrametricEvaluation.r, ...



What are its inputs. 

What its outputs. 

Which dependencies is has. 

**Which scripts it runs automatically.**
The following r libraries:
- SimSpeciesNetworks.r
- Net_HybridLamb_Format.jl
- UltrametricEvaluation.r

## `SimSpeciesNetworks.r`

What the script does. 

What are its inputs. 

What its outputs. 

**Which dependencies is has.**
- SiPhyNetwork
- geiger

Which scripts it runs automatically.
nothing

## `Net_HybridLamb_Format.jl`

What the script does. 

What are its inputs. 

What its outputs. 

Which dependencies is has. 

Which scripts it runs automatically.

## `UltrametricEvaluation.r`

What the script does. 

What are its inputs. 

What its outputs. 

Which dependencies is has. 

Which scripts it runs automatically.

## `SnaqHybridzation.jl`

What the script does. 

What are its inputs. 

What its outputs. 

Which dependencies is has. 

Which scripts it runs automatically.
