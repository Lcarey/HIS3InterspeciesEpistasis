README file for evfold.org results
----------------------------------------------------------------------------

Files are organized into subdirectories. Depending on the job configuration,
subdirectories may include:

alignment : alignment files (constructed, supplied, or retrieved)
contact_maps : constraint plots and contact map plots
ev_couplings : residue pair scoring and ranking / selection / summary
job_config : settings chosen for this job
residue_numbering : mapping between uniprot, alignment, pdb positions
sequence : complete proteins sequence and secondary structure analysis
models_compared_to_known_structure : comparisons to known structure from pdb
structure_inputs : ECs and secondary structure inputs to cns_solve
structure_outputs : computed structures and non-comparitive rankings

Further Details of directory contents:

Subdirectory: alignment
----------------------------------------------------------------------------
File: {jobname}_{domainname}_user_supplied_alignment.fa
    the uploaded alignment, when present
File: {jobname}_{domainname}.fa
    the fasta format alignment file used to calculate ECs
File: {jobname}_{domainname}.weight_table
    the number of closely related sequences for each alignment member and
    corresponding re-weighting score
File: {jobname}_{domainname}_{tool}_e_value_table.csv
    alignment statistics for alignment generation runs at different e-values,
    including number of sequences, effective number of sequences,
    domain width, counts of "match state" and "non-match state" columns,
    coverage, and scores for different evalues

Subdirectory: contact_maps
----------------------------------------------------------------------------
File: {jobname}_ConstraintMap_{X}.fig
    matlab figure file containing a constraint map showing X ECs -- a 2D
    representation of evolutionary couplings between residues within the
    protein (orange-red asterisks represent ECs)
File: {jobname}_ConstraintMap_{X}.pdf
    as previous, but in PDF format
File: {jobname}_ContactMap_{X}.fig
    matlab figure file containing a contact map showing X ECs -- a 2D
    representation of evolutionary couplings between residues within the
    protein (orange-red asterisks represent ECs, grey represent known pdb
    structure contacts, light blue represent regions not in the known
    pdb structure, green represent computed structure contacts)
File: {jobname}_ContactMap_{X}.pdf
    as previous, but in PDF format

Subdirectory: ev_couplings
----------------------------------------------------------------------------
File: {jobname}_{scoringmethod}.txt
    list of all-by-all residue pairings, and score computed by chosen method
    MI_DI column headers:
    - 1stResidueNum
    - 1stResidueCode
    - 2ndResidueNum
    - 2ndResidueCode
    - mutual information score
    - DI score
    PLM columns are the same, replacing DI score with PLM score, and
    omitting MI scores (always 0)
File: {jobname}_CouplingScores.csv
    sorted (descending) EC scores with filtration codes. Columns headers:
    - 1stResideNum
    - 2ndResidueNum
    - ECscore
    - placeholder (always 0)
    - 1stResidueAlignColumnConservation
    - 2ndResidueAlignColumnConservation
    - secondary structure fiter (flag "999")
    - high conservation filter (flag "888")
    - cys-cys filter (flag "222")
    - residue code for 1st Residue
    - residue code for 2nd residue
    Note: when a comparison known pdb structure is given, the placeholder
    column is filled with a euclidian distance from the known structure,
    and the resulting file will be located in the
    models_compared_to_known_structure subdirectory
File: {jobname}_EC_residues_summary_table.csv
    EC scores summed for each residue, sorted in descending order. Column
    headers are:
    - residue index
    - amino acid
    - number of ECs
    - cumulative strength
    - EC strength
    - conservation
File: {jobname}_FP_Plot.fig
    A "false positive plot" showing the euclidian distance in the known
    pdb structure against the EC score rank. ECs at high distances might be
    regarded as false positives
File: {jobname}_FP_Plot.pdf"
    as previous, but in PDF format
*See note below for additional visualization files

Subdirectory: job_config
----------------------------------------------------------------------------
File: {jobname}_user_config.txt":
    user supplied configuration settings

Subdirectory: residue_numbering
----------------------------------------------------------------------------
File: {jobname}_{protein_id}.indextable
    secondary structure predictions and multiple sequence alignment
    conservation.
    Column headers are:
    - uniprot index
    - uniprot residue code
    - secondary structure prediction
    - secondary structure confidence
    - alignment column index
    - alignment column conservation pct
    - alignment columm match state flag: * = match column, ~ = non-match
    - in construct flag (residue would be part of cns_solve modeling)
File: {jobname}_{protein_id}.indextableplus
    contains all the information in the indextable, with these additional columns
    pertaining to PDB mappings:
    - pdb_atom
    - pdb_chain
    - pdb_index
    - pdb_residue
    - pdb_x_pos
    - pdb_y_pos
    - pdb_z_pos

Subdirectory: sequence
----------------------------------------------------------------------------
File: {jobname}_{protein_id}.fa *or* {jobname}_user_supplid_sequence.fa
    the protein sequence of interest in FASTA format
File: {jobname}_{protein_id}_domain_region.fa
    the sequence of interest, trimmed to correspond to the domain region
    specified, in FASTA format
File: {jobname}_psipred.txt
    psipred secondary structure predictions. column headers:
    - residue number
    - residue code
    - secondary structure feature assignment code E=strand H=helix C=coil
    - secondary structure feature assignment confidence
File: {jobname}_topology3state.txt
    memsat-svm (or user provided) secondary structure predictions.
    column headers:
    - residue number
    - residue code
    - secondary structure feature assignment code I=inward helix edge
        X=helix middle, O=outward helix edge (for TM helices)
    - secondary structure feature assignment confidence

Subdirectory: models_compared_to_known_structure
----------------------------------------------------------------------------
File: {jobname}_pdb{XXXX}.ent
    the reference known structure from PDB for comparison
File: {jobname}_CouplingScoresCompared.csv
    identical to ev_couplings/{jobname}_CouplingScores.csv except that the
    placeholder column is filled in with the C-alpha/NearestAtom Euclidian
    Distance between residues
File: {jobname}_{N}_pymol.txt
    rmsd calculations by PyMol (one file per bin, N = bin constraint count)
File: {jobname}_pymol_comparison_summary.csv
    summary of PyMol RMSD results. Column headers are:
    - computed_structure
    - RMSD
    - matched_atoms (residues)
    - alignment_cycles

Subdirectory: structure_inputs
----------------------------------------------------------------------------
File: {jobname}_{N}_Couplings.tbl
    selected EC constraints input for CNS (one file per bin, N = count)
File: {jobname}_SS_distance.tbl
    secondary structure distance constraints input for CNS
File: {jobname}_SS_angle.tbl:
    secondary structure torsion angle constraint input for CNS
File: {jobname}.seq
    the domain modeling sequence (from the sequence of interest), formatted
    for input to CNS

Subdirectory: structure_outputs
----------------------------------------------------------------------------
File: {jobname}_{N}_{V}.hMIN.pdb
    computed structure files after energy minimization. N = count of applied
    EC constraints, V = structure variant number
File: {jobname}_alphabeta_ranking.txt or {jobname}_transmembrane_ranking.txt
    the (non-comparitive) ranking file, with a score for each structure
    based on inherent properties and extent of constraint satisfaction
    Column headers for alphabeta_ranking:
    - filename
    - #constraints
    - variant
    - #alpha_residues
    - #alpha_dihedrals
    - #good_alpha_dihedrals
    - good_alpha_dihedral_fraction
    - #beta_residues
    - #beta_dihedrals
    - good_beta_proportion
    - total_weighted_score
    Column headers for transmembrane_ranking:
    - filename
    - bin (number of EC constraints applied)
    - num (variant number)
    - lipid (hydrophobic exposure)
    - ss-agr (secondary structure agreement (sequence / structure))
    - em (satisfaction of EC constraints)
    - total
File: {jobname}_rank{R}_{N}_{V}.hMIN.pdb
    These files are identical to the matching computed structure files,
    but with the rank index R added to the filename to make sorting
    and selecting top ranking computed structures easier.
File: {jobname}_draw_ss.pml
    a script for drawing the secondary structure on a computed structure
*See note below for additional visualization files

*note: The following three results scripts will either be located in the
    ev_couplings subdirectory (for EVfold/structure computation jobs) or
    the structure_outputs subdirectory (for EVcouplings/non-structure
    computation jobs). When located in the structure_outputs subdirectory,
    they visualize ECs on the computed structures (which must first be loaded
    into pymol before running the visualization script). Otherwise, they
    visualize ECs on the known PDB structure, which will be loaded
    automatically by the script.
File: {jobname}_draw_EC_lines.pml
    a PyMol script for drawing ECs as lines between resides
File: {jobname}draw_EC_residues_sausage.pml
    a PyMol script for drawing EC cumulative strength (thicker sections
    indicate higher EC cumulative strength)
File: {jobname}_draw_EC_residues_spheres.pml
    a PyMol script for drawing resides with strong cumulative EC strength
    as colored spheres (Red: strongest, Orange: strong)

----------------------------------------------------------------------------
