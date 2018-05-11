%% Reviews for Nature Ecology & Evolution. May 2018
% LBC

% Reviewer #1 (Remarks to the Author):

%When reviewing the original submission, my main concern was that fitness
%estimates may be inaccurate due to the low frequency of individual
%variants, the high similarity between variants, and a potentially high
%frequency of sequencing errors. A typical variant is present in its
%sequencing library at a frequency of 1:300,000, and it is critical that
%the authors show experimental results to prove their ability to accurately
%quantify these rare variants. The authors have provided theoretical
%arguments and experimental controls that support the accuracy of their
%measurements, but these arguments do not entirely dispel my doubts, and
%raise additional questions.

% -- prove their ability to accurately quantify these rare variants

%The authors argue that the probability that a read contains a sequencing
%error is the important parameter to estimate, and the relative frequency
%of variants and sequencing errors is secondary. I disagree. Imagine, for
%example, an experiment in which the entire library consisted of one
%variant: wild-type HIS3. In this experiment, even if 90% of the reads were
%accurate, the remaining 10% of reads would contain sequencing errors, that
%the authors would identify as mutated variants of HIS3. Clustering of
%these erroneous reads using Starcode would only partially mitigate the
%problem. Importantly, even though 90% of all reads were accurate, 100% of
%the mutated variants identified would be false, which shows that there is
%no simple relationship between the proportion of error-free reads and the
%proportion of error-free variants identified.

%Further, the authors use a binomial model to conclude that only 10% of
%reads contain errors. This calculation relies on the assumption of 0.1%
%per-nucleotide sequencing error rate, and of equal distribution of
%sequencing errors across reads. Both assumptions may be incorrect. For
%example, the Starcode paper states that Illumina has a 1–2% error rate
%consisting of substitutions. Without careful controls it is impossible to
%tell what is the error rate in any particular experiment, and the
%frequency of erroneous reads may be higher than 10%.

%The authors also present a spike-in experiment. The spike-in figures
%appear to show that each spike-in variant was present at a frequency
%between 0.05-0.3 (5%-30%), much higher than sequencing error rates. It is
%good to see that such frequent variants can be quantified with reasonable
%accuracy, but it is still unclear whether the same applies to rare, highly
%similar variants. Thus, the spike-in experiment does not address my main
%concern. The authors also mention a number of internal controls that
%support the vailidity of their experiment, including the fitness estimates
%of frameshift and synonymous variants. I agree that these controls support
%several of the conclusions, but they do not really address my main
%concern, which is the accuracy of fitness estimation of rare variants.

% --  the accuracy of fitness estimation of rare variants.

%The authors also state that “Starcode converts the variants at distance 1
%or 2 from a more abundant variant to this more abundant variant”, which
%raises additional questions. First, this filtering makes it difficult to
%measure the fitness of variants that are one amino acid substitution away
%from each other, because many such variants are separated by 1 or 2
%nucleotide substitutions. The authors report the measurement of many
%variants that are one amino acid substitution away from each other, for
%example in Figures 2a, b. I do not understand how this analysis is
%compatible with the Starcode filtering step mentioned above. Second, the
%Starcode filtering step will artificially increase the read counts
%associated with certain variants, at the expense of other variants.
%Without careful analysis of this effect, it is difficult to tell how this
%may influence the conclusions, for example regarding epistasis, or
%regarding the dependence of the effect of individual substitutions on the
%genetic background.

% -- measure the fitness of variants that are substitution apart vs Starcode clustering

%% LBC
% use synonymous variants to measure our ability to correctly call fitness
% for low frequency genotypes

%% load synonymous data data
T = readtable('~/Develop/HIS3InterspeciesEpistasis/Data/synonymous_variants_rescaled_data.tab',...
    'Delimiter','\t','Filetype','text','format','%d%s%s%f%f');


%% %%  %% %%  %% %%  %% %% %% Part 1 %% %%  %% %%  %% %%  %% %%
% -- prove their ability to accurately quantify these rare variants
% --  the accuracy of fitness estimation of rare variants.

%% load data and calculate the distance between fitness measurements for synonymous variants and mean across all syn variants
s = struct();
for SegN = 1:12
    SegN
    save( 'freq_effect.mat' , 's');
    s(SegN).SegN = SegN ;
    Q = T(T.SegN == SegN ,:);
    G = grpstats( Q , 'aa_seq' , {'mean' 'median' 'std' 'sem'} , 'DataVars' , {'t0_fr' 's'});
    G2 = G(G.GroupCount>10,:);
    
    data = NaN( 1e6 , 3);
    c = 1 ;
    for I = 1:height(G2)
        idx =  ismember(Q.aa_seq,G2.aa_seq{I});
        syn_fit = Q.s(idx);
        syn_frq = Q.t0_fr(idx);
        n = sum(idx)-1;
        idx2 = c:(c+n);
        data(idx2,1) = syn_frq ;
        data(idx2,2) = syn_fit ;
        data(idx2,3) = repmat(G2.mean_s(I),n+1,1) ;
        c = c+n ;
    end
    data = data(1:c,:) ;
    s(SegN).data = data ;
    %
    Y = abs(data(:,2)-data(:,3)) ./ 0.45 ;
    X = data(:,1) ./ sum(data(:,1)) ;
    XGRP = zeros(size(X));
    XGRP( X > 1/500000 & X < 1/100000) = 1;
    XGRP( X > 1/100000 & X < 50000) = 2;
    XGRP( X > 1/50000  & X < 10000) = 3;
    XGRP( X > 1/10000 ) = 4;
    s(SegN).X = X;
    s(SegN).XGRP = XGRP;
    s(SegN).Y = Y;
end
%% ecdfs for each segment
lbls = { '<1/500,000' '<1/100,000' '<1/50,000' '<1/10,000' '>1/10,000'} ;
for SegN = 1:12
    try
        ux = unique(s(SegN).X);
        fh = figure('units','centimeters','position',[5 5 8  8]);
        hold on ;
        for I = 1:numel(ux)
            [f,x] = ecdf( s(SegN).Y( s(SegN).X==ux(I)));
            plot(x,f,'LineWidth',3);
        end
        ylabel('Fraction of genotypes')
        xlabel('Difference in fitness')
        grid on ;
        legend( lbls(ux+1) , 'location','best')
        set(gca,'ytick',0:.1:1)
        title( ['Segment ' num2str(SegN)]);
        xlim([0 0.75])
    catch
    end
end

%% boxplot of above data
data = NaN(0);
for SegN = 1:12
    data = vertcat( data , [ s(SegN).X s(SegN).Y] );
end
fh = figure('units','centimeters','position',[5 5 14 6]);
boxplot(data(:,2) , data(:,1) ,'symbol','')
ylim([0 0.5])
set(gca,'xticklabel',lbls)
ylabel('Difference in fitness')
set(gca,'ytick',0:.1:1)


%% %%  %% %%  %% %%  %% %% %% Part 1b %% %%  %% %%  %% %%  %% %%
% -- prove their ability to accurately quantify these rare variants
% --  the accuracy of fitness estimation of rare variants.
%    for nonsense mutants
%% Fitness of nonsense variants, grouped by t0_fr
idx = regexpcmp(T.aa_seq,'_');
X = T.t0_fr ./ 1e6 ;
XGRP = zeros( height(T),1);
XGRP( X > 1/500000 & X < 1/100000) = 1;
XGRP( X > 1/100000 & X < 50000) = 2;
XGRP( X > 1/50000  & X < 10000) = 3;
XGRP( X > 1/10000 ) = 4;
%%
fh = figure('units','centimeters','position',[5 5 14 6]);
boxplot(T.s(idx) , XGRP(idx) ,'symbol','')
ylim([0 0.5])
set(gca,'xticklabel',lbls)
ylabel('Difference in fitness')
set(gca,'ytick',0:.1:1)

ug = unique(XGRP);
fh = figure('units','centimeters','position',[5 5 8  8]);
hold on ;
for I = 1:4
    [f,x] = ecdf( T.s(idx & XGRP==ug(I)));
    plot(x,f,'LineWidth',3);
end
ylabel('Fraction of genotypes')
xlabel('Difference in fitness')
grid on ;
legend( lbls(1:4) , 'location','best')
set(gca,'ytick',0:.1:1)
xlim([0 0.5])


%% %%  %% %%  %% %%  %% %% %% Part 2 %% %%  %% %%  %% %%  %% %%
% % -- measure the fitness of variants that are substitution apart vs Starcode clustering
s2 = struct();
for SegN = 1:12
    SegN
    sT = readtable( ['~/Develop/HIS3InterspeciesEpistasis/Data/S' num2str(SegN) '_scaled_info_v2.csv'] , 'Delimiter','\t');
    sT = sT( sT.len == mode(sT.len) ,:);
    
    nonsense_mutants = sT.aa_seq( logical(sT.nonsense)) ;
    nonsense_mutants_fit = sT.s( logical(sT.nonsense)) ;
    fit_mutants = sT.aa_seq(sT.size>=3 & sT.s>0.7) ;
    fitness_1 = NaN(0);
    fitness_2 = NaN(0);
    fitness_3 = NaN(0);
    fitness_4 = NaN(0);
    fitness_5 = NaN(0);
    
    
    for I = 1:numel(fit_mutants)
        d = cellfun( @(X)HammingDistance( X ,fit_mutants{I}) , nonsense_mutants);
        fitness_1 = vertcat( fitness_1 , nonsense_mutants_fit(d==1));
        fitness_2 = vertcat( fitness_2 , nonsense_mutants_fit(d==2));
        fitness_3 = vertcat( fitness_3 , nonsense_mutants_fit(d==3));
        fitness_4 = vertcat( fitness_4 , nonsense_mutants_fit(d==4));
        fitness_5 = vertcat( fitness_5 , nonsense_mutants_fit(d==4));
        
    end
    
    s2(SegN).nonsense_mutants = nonsense_mutants ;
    s2(SegN).nonsense_mutants_fit = nonsense_mutants_fit ;
    s2(SegN).fit_mutants = fit_mutants ;
    s2(SegN).fitness_1 = fitness_1 ;
    s2(SegN).fitness_2 = fitness_2 ;
    s2(SegN).fitness_3 = fitness_3 ;
    s2(SegN).fitness_4 = fitness_4 ;
    s2(SegN).fitness_5 = fitness_5 ;
    
end
%%
fh = figure('units','centimeters','position',[5 5 8  8]);
hold on ;
for I = 1:4
    [f,x] = ecdf( eval(['fitness_' num2str(I)]) );
    plot(x,f,'LineWidth',3,'DisplayName',num2str(I));
end
ylabel('Fraction of genotypes')
xlabel('Fitness of nonsense mutants')
grid on ;
legend( 'location','best')
%set(gca,'ytick',0:.1:1)
%xlim([0 0.5])
