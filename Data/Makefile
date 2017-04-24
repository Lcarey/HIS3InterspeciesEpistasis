# process the data files to contain only aa-seq, lib-type and fitness
# LBC 2017
smaller_data_files:
	for N in 1 2 3 4 5 6 7 8 9 10 11 12; do \
		cut.pl -n aa_seq,s,len,nogap,stop,nonsense,middle,nat_lib,nat,lib,size < S$${N}_scaled_info.csv > $${N}_all.tab ;\
		cut.pl -n aa_seq,s,len,nogap,stop,nonsense,middle,nat_lib,nat,lib,size < S$${N}_scaled_info.csv | grep '\t1\t0\t0\t1\t1\t1\t1\t' > $${N}_natlib.tab ;\
		cut.pl -n aa_seq,s,len,nogap,stop,nonsense,middle,nat_lib,nat,lib,size < S$${N}_scaled_info.csv | grep '\t1\t0\t0\t1\t1\t[01]\t[01]\t' > $${N}_nogap_nostop_middle.tab ;\
	done

syndata:
	rm -f syndata.tab x
	for N in 1 2 3 4 5 6 7 8 9 10 11 12; do \
		cut.pl -n aa_seq,s < Synonymous/S$${N}_sum_10_size_4_fit.csv | add_column.pl -b -s $${N} >> x ;\
	done
	grep -v aa_seq x > syndata.tab
	rm -f x



