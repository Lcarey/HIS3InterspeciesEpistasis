# process the data files to contain only aa-seq, lib-type and fitness
# LBC 2017
smaller_data_files:
	for N in 1 2 3 4 5 6 7 8 9 10 11 12; do \
		cut.pl -n aa_seq,s,len,nogap,stop,nonsense,middle,nat_lib,nat,lib,size < S$${N}_scaled_info.csv > $$N.tab ;\
	done

