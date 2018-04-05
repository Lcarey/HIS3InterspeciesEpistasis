%% We performed least square fit to an exponential growth model on original frequencies, rather than to linear model on log_frequencies because exponential growth model could correctly deal with genotypes that drop out of the competition. X% of genotypes and X% of extant genotypes have fr_t0>0, fr_t1>0, fr_t2=0. Exponential growth model allows us analyse them throw the same procedure as genotypes with ft_t0,fr_t1,fr_t2>0, while if work with log ratios there are no good ways to combine information of genotype fitness in t0-t1 interval and information of the fact that it drops out of the competition in t1-t2 interval.
% LBC April 2018

%calc for each segment
DATADIR = '~/Develop/HIS3InterspeciesEpistasis/Data/' ;
R = dataset() ;
R.SegN = NaN(12,1) ;
R.N_strains_all_genotypes  = NaN( 12,1) ;
R.N_strains_dropout_all_genotypes = NaN( 12,1);
R.N_strains_extant_genotypes  = NaN( 12,1) ;
R.N_strains_dropout_extant_genotypes = NaN( 12,1);
 
R.N_strains_nonsense_genotypes  = NaN( 12,1) ;
R.N_strains_dropout_nonsense_genotypes = NaN( 12,1);


for SegN = 1:12
    T = readtable( [DATADIR 'S' num2str(SegN) '_scaled_info_v2.csv'] , 'FileType','text','Delimiter','\t');
    idx_dropout_2 = T.t2_fr==0 & T.t1_fr > 0 ;
    idx_extant = T.nat_lib & ~T.stop & T.nogap & ~T.nonsense ;
    idx_nonsense =   T.stop & T.nonsense ;
    
    R.N_strains_all_genotypes(SegN) = height(T) ; 
    R.N_strains_dropout_all_genotypes(SegN) = sum(idx_dropout_2) ;
    
    R.N_strains_dropout_extant_genotypes(SegN) = sum(idx_dropout_2 & idx_extant) ;
    R.N_strains_extant_genotypes(SegN) = sum(idx_extant) ;

    R.N_strains_dropout_nonsense_genotypes(SegN) = sum(idx_dropout_2 & idx_nonsense) ;
    R.N_strains_nonsense_genotypes(SegN) = sum(idx_nonsense) ;

    
    R.SegN(SegN) = SegN ; 
    
    R.MeanFit__Dropout_ExtantGenotypes(SegN) = mean(T.s(idx_dropout_2 & idx_extant)) ; 
    R.MeanFit__Dropout_AllGenotypes(SegN) = mean(T.s(idx_dropout_2 )) ; 
    R.MeanFit__NOTDropout_ExtantGenotypes(SegN) = mean(T.s(~idx_dropout_2 & idx_extant)) ; 
    R.MeanFit__NOTDropout_AllGenotypes(SegN) = mean(T.s(~idx_dropout_2 )) ; 
end

%%
R.pct_dropout_all_genotypes = R.N_strains_dropout_all_genotypes ./ R.N_strains_all_genotypes * 100 ; 
R.pct_dropout_extant_genotypes = R.N_strains_dropout_extant_genotypes ./ R.N_strains_extant_genotypes * 100 ; 
R.pct_dropout_nonsense_genotypes = R.N_strains_dropout_nonsense_genotypes ./ R.N_strains_nonsense_genotypes * 100 ; 

writetable( dataset2table(R)  , '~/Downloads/dropout_t2_statistics.xlsx');

%%
sum(R.N_strains_dropout_all_genotypes) ./  sum(R.N_strains_all_genotypes)  * 100 

sum(R.N_strains_dropout_extant_genotypes) ./  sum(R.N_strains_extant_genotypes)  * 100 

sum(R.N_strains_dropout_nonsense_genotypes) ./  sum(R.N_strains_nonsense_genotypes)  * 100 