% Fig_CalcDiversityAcrossMultipleSequenceAlignment
% To study the evolutionarily-relevant region of the sequence landscape we selected His3, a gene coding for imidazoleglycerol-phosphate dehydratase (IGPD, His3p), a conditionally-essential enzyme and identified XXX amino acid differences between His3 orthologues from 21 different yeast species spanning 500 million years of sequence divergence (Rhind 2011; Figure 1B and 1C). XXX amino acid substitutions correspond to XXX possible genotypes, too many to survey with modern methods.
%
% May 11, 2017 LBC

DDIR = '~/Develop/HIS3InterspeciesEpistasis/Data_Small_Tables/';
%BigMSA = fastaread([DDIR 'HIS3_569species_uniprot_alignment.fas']);
BigMSA = fastaread([DDIR 'aa_seq.txt']);

mat = vertcat(BigMSA.Sequence);
n_possibilities = NaN( length(mat(1,:)) ,1);
for I = 1:length(n_possibilities)
    aas_at_this_pos =  unique(   mat(:,I)  ) ;
    n_possibilities(I) = sum(~ismember( aas_at_this_pos ,'-')) ;
end

cp = cumprod(n_possibilities);
cs = sum(n_possibilities(n_possibilities>1))
cp(end)
%%
figure;
plot( n_possibilities,'-k')
axis tight;
ylabel('# different AAs')
xlabel('Position in multiple sequence alignment')
title( sprintf( '%0.05e possible AA seqs' , cp(end)))