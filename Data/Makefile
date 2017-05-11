# process the data files to contain only aa-seq, lib-type and fitness
# LBC 2017
smaller_data_files:
	for N in 1 2 3 4 5 6 7 8 9 10 11 12; do \
		cut.pl -n aa_seq,s,len,nogap,stop,nonsense,middle,nat_lib,nat,lib,size < S$${N}_scaled_info_v2.csv > $${N}_all.tab ;\
	done

syndata:
	rm -f syndata.tab x
	for N in 1 2 3 4 5 6 7 8 9 10 11 12; do \
		cut.pl -n seq,aa_seq,s,t0_fr, < Synonymous/S$${N}_nucl_scaled.csv | add_column.pl -b -s $${N} >> x ;\
	done
	echo 'SegN\tseq\taa_seq\ts\tt0_fr' > syndata.tab
	grep -v aa_seq x >> syndata.tab
	rm -f x
	mv syndata.tab synonymous_variants_rescaled_data.tab



