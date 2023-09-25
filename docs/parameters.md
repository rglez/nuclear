# Parameters

NUCLEAR uses a configuration file to specify all their required parameters. Users must decide what job to perform:
(i) a molecular hotspot determination or (ii) an oligonucleotide search.

For a **hotspot search**, the sections [inputs], [mcss_files], and [spots_params] are the only ones that should 
appear in the configuration file. If an **oligonucleotide search** is requested, then the mandatory sections are
[inputs], [mcss_files], [zones_params], [search_params], and [minimization].

Below, all section parameters have been described.



## **Section \[inputs\]**

- **input_dir:** Path to the directory containing MCSS explorations   
 *Acceptable Values:* A valid path  
 *Recommended Value:* None  
 *Description:* Specifies the path to find the *.crd* files of MCSS docking results. NUCLEAR will check if this path exists.

  
- **output_dir:** Path to the output directory   
  *Acceptable Values:* A valid path 
  *Recommended Value:* None  
  *Description:* Specifies the root path where NUCLEAR writes results. NUCLEAR will check if this path exists, in which case the program will stop to avoid overwriting previous results.  
  

- **prot_coords:** Path to the *receptor.pdb*  
  *Acceptable Values:* A valid path  
  *Recommended Value:* None  
  *Description:* Specifies the path to the *receptor.pdb* coordinates file. NUCLEAR will check if this path exists and that the receptor file format is *.pdb*.  
  

- **ncores:** Number of cores for parallelization   
  *Acceptable Values:* integer >= 1  
  *Recommended Value:* N (where N is the maximum number of cores available)  
  *Description:* Specifies the number of cores used to construct the matrix of clashes.  
  


## **Section \[mcss_files\]**   
  
\[mcss_files\] is a mandatory section that contains no parameter. Instead, users must specify one or more tuples of three blank-separated columns. Note that one of the two last columns (or both) must be N (unrestricted) to avoid ambiguities in the selection of the fragments' poses from docking distributions.  
  
**1st column (string):** The name of the *.crd* file (extension included). Only include the name, but not the path, which was already specified in the `input_dir` parameter.
  
**2nd column (float):** Maximum score to be considered in the corresponding MCSS distribution. Fragments' poses having a higher score will not be considered.
  
**3rd column (integer):** Number of poses to be considered in the corresponding MCSS distribution. 
  


## **Section \[spots_params\]**  
  
- **prot_sel** NUCLEAR receptor selection   
  *Acceptable Values:* \[ all | heavy | protein \]  
  *Recommended Value:* protein  
  *Description:* Selection of receptor atoms considered when calculating receptor-ligand contacts. Options are limited to `all`: all atoms in the `prot_coords` file, `heavy`: non-hydrogen atoms in the `prot_coords` file, and `protein`: protein atoms in the `prot_coords` file.  
  

- **frag_sel** NUCLEAR fragment (ligand) selection   
  *Acceptable Values:* \[ all, all_nopatch, heavy, heavy_nopatch \]  
  *Recommended Value:* all_nopatch  
  *Description:* Selection of fragment atoms considered when calculating receptor-ligand contacts. Options are limited to `all`: all atoms in the *ligand.crd* file, `all_nopatch`: all atoms except patch ones in the *ligand.crd* file, `heavy`: non-hydrogen atoms in the *ligand.crd* file, and `heavy_nopatch`: non-hydrogen atoms except patch ones in the *ligand.crd* file. Patch atoms (not used when constructing oligonucleotides) are defined as:*\['C5T', 'O5T', 'O1P', 'O2P', 'H5T1', 'H5T2', 'H5T3'\]*.  
   

- **inter_dist** Inter-atomic distance (in Angstroms) for contacts calculation   
  *Acceptable Values:* float > 0.0  
  *Recommended Value:* 3.5  
  *Description:* Inter-atomic distance under which a contact between a ligand atom and a protein atom is declared.  
  

- **density_cut** In hotspots detection, clusters of fragments' poses having a density value under this cutoff get filtered out   
  *Acceptable Values:* 0.0 < float < 1.0  
  *Recommended Value:* 0.05  
  *Description:* Clusters having a density value under this cutoff get filtered out. The density of a cluster is calculated as the number of poses it contains divided by the number of distinct residues involved in the contact. This magnitude is then normalized to the 0-1 interval.  
  

- **merge_cut** In hotspot detection, cluster seeds having a similarity value under this cutoff provokes their clusters to merge   
  *Acceptable Values:* 0.0 < float < 1.0  
  *Recommended Value:* 0.25  
  *Description:* In hotspot detection, cluster seeds having a similarity value under this cutoff provokes their clusters to merge. The similarity metric employed is the Tanimoto index (lowest values correspond to more similar seeds).  
   



## **Section \[zones_params\]**  
  
- **prot_sel**< idem to what is described in section **\[spots_params\]**>  
- **frag_sel**< idem to what is described in section **\[spots_params\]**>  
- **inter_dist**< idem to what is described in section **\[spots_params\]**>  
  


## **Section \[search_params\]**  
  
- **top_N** Number of top-N best-scored sequences to be returned   
  *Acceptable Values:* integer >= 1  
  *Recommended Value:* 1000  
  *Description:* Number of top-N best-scored sequences to be returned. The sequence score is calculated as the summation of the scores of each composing nucleotide.  
   

- **seq_min_size** Minimum number of nucleotides in returned sequences   
  *Acceptable Values:* integer >= 2  
  *Recommended Value:* 2  
  *Description:* Minimum number of nucleotides in returned sequences. When `seq_to_search` is not `all`, this value must equal the number of nucleotides specified in the `seq_to_search` parameter.  
  
  - **seq_max_size** Maximum number of nucleotides in returned sequences   
  *Acceptable Values:* integer >= `seq_min_size`  
  *Recommended Value:* N  
  *Description:* Maximum number of nucleotides in returned sequences. Users can retrieve the maximum number of linkable nucleotides by setting this to N. When `seq_to_search` is not `all`, this value must equal the number of nucleotides specified in the `seq_to_search` parameter.  
  

- **seq_to_search** Sequence of nucleotides to search   
  *Acceptable Values:* \[ all, colon-separated string of resnames \]  
  *Recommended Value:* all  
  *Description:* Sequence of nucleotides to search. If not `all`, the oligonucleotide sequence specified must comprise colon-separated strings corresponding to each nucleotide resname, as found in the corresponding *.crd* files. The underscore _ specifies any nucleotide. When `seq_to_search` is not `all`, parameters `seq_max_size` must equal `seq_min_size`, and both must be set to the Number of nucleotides specified in `seq_to_search`.  
  

- **path_to_search** Atomic selection of receptor's residues to search   
  *Acceptable Values:* valid string selection (defined in agreement to ProDy syntax)  
  *Recommended Value:* all
  *Description:* Atomic selection of receptor residues. Only contacts containing these residues indices will be considered. 
   

- **max_dist_O3C5** The maximum distance between O3' and C5' atoms that may link consecutive nucleotides.   
  *Acceptable Values:* float > 0.0  
  *Recommended Value:* 4.5 A  
  *Description:* The maximum distance between O3' and C5' atoms that may link consecutive nucleotides.  
   

- **clash_dist** Minimum distance to declare an atomic clash   
  *Acceptable Values:* 0.0 < float < `max_dist_O3C5` 
  *Recommended Value:* 2.0 A  
  *Description:* Minimum distance to declare an atomic clash. When atoms of different nucleotides inside an under-construction sequence are closer than this distance, a clash is declared, and the nucleotide provoking the clash is discarded from the sequence.  
   


## **Section \[minimization\]**  

- **minimize** Produce CHARMM input files for minimization?  
  *Acceptable Values:* [True, False]
  *Recommended Value:* True 
  *Description:* If set to True, this option makes NUCLEAR to produce CHARMM's minimization input for every oligonucleotide sequence produced.   
   

- **prot_topol** Path to the CHARMM topology file required for minimization  
  *Acceptable Values:* A valid path
  *Recommended Value:* None 
  *Description:* Path to the CHARMM topology file required for minimization.   
   

- **prot_param** Path to the CHARMM parameter file required for minimization  
  *Acceptable Values:* A valid path
  *Recommended Value:* None 
  *Description:* Path to the CHARMM parameter file required for minimization.   
