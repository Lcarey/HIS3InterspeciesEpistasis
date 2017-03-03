aa_seq	Sequence of genotype	
		It consists of first mutated part, middle not-mutated part, and last mutated part.
		These parts are ~10 amino acids in length but not exactly. 
		Look in table positions.csv for more information about design of mutated parts

tX_fr	Frequencies of genotype in timepoint tX, 10^-6

tX_fr_var	Variance of genotype frequency 

size	Number of nucleotide genotypes corresponding to amino acid sequence

s		Fitness 
		Estimated as exponent fitted in growth curve.
		Fequencies of nucleotide genotypes corresponding to one amino acid genotype were summed before fitting. 
		Scaled between 0 and 0.45, where 0 is fitness of nonsense genotypes summed together and 0.45 is wt fitness

s_std	Error of fitness reported by fitting function.

y0		Initial frequency of genotype, *10^-6
		Estimated as a coefficient when fitting exponent in growth curve
		
y0_std	Error of y0 reported by fitting function

s_original - fitness before scaling

s_std_original - fitness eror before scaling
		
len		Length of aa_seq

nonsense=1 means that there is either frame shift or stop codon

nogap	=1 means that there are no gaps when align aa_seq of reference sequence

middle	=1 means that there are no mutations in middle part which is by design mutation free

nat		=1 means that aa_seq consists only of mutations which appear in other species (based on Vika's alignment of 569 species)

lib		=1 means that aa_Seq consists only of mutations which were in designed library

nat_lib	=1 means that (lib==1)&(nat==1)

mut_list	List of all mutations in aa_seq relative to Scer	
			Format 0G:5T:23A
			Where 0,5,23 are positions of mutation reltive to segment (starts from 0) and G,T,A are new amino acids		
			It was calculated only for (nopgap==1)&(middle==1)
			
dist_X	Hamming distance between genotype and corresponding species

dist_node_X		Hamming distance between genotype and internal node of philogenetic tree 

dist_min_all	Minimal among all distances

dist_min_sp		Minimal among distances to species (not internal tree nodes)
			
