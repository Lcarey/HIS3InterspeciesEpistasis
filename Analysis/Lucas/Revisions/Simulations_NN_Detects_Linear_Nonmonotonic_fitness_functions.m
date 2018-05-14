%% get fitness potential values for all AA states
ltrs = 'AGLMFWKQESPVICYHRNDT';
ltrs_fit_pot = random( 'normal' , 0 , 0.5 , length(ltrs), 1);


fh = figure('units','centimeters','position',[5 5 8 8]);
[~,o]=sort(ltrs_fit_pot,'descend');
plot( 1:20 , ltrs_fit_pot(o) , '.w')
text( 1:20 , ltrs_fit_pot(o) , ltrs(o)')
set(gca,'xtick',[])
ylabel('Effect on fitness potential')
xlim([0 21])

T=table();
T.AA = ltrs' ;
T.FP = ltrs_fit_pot ; 
writetable( T , '~/Develop/HIS3InterspeciesEpistasis/Analysis/Lucas/Revisions/Simulations_NN_Detects_Linear_Nonmonotonic_fitness_functions__AAstates.tab' ...
    , 'FileType','text','Delimiter','\t');
%% build up genotypes and fitness potential
genotype_length = 20 ; 
Nseqs = 5e5 ; 
max_FP = 6 ;
min_FP = -3 ; 

aa_seq = cell(  Nseqs , 1 );
fitness_potential  = NaN( Nseqs , 1 );

for I = 1:Nseqs
	idx = randi(20 , genotype_length , 1); 
	aa_seq{I} = ltrs(idx) ;
	fitness_potential(I) = sum( ltrs_fit_pot(idx) ) ;
end


% we sparely populate high FP values. get some more
aa_seq2 = cell( Nseqs , 1);
fitness_potential2 = NaN( Nseqs , 1);
for I = 1:Nseqs
	idx = randi(20 , genotype_length , 1); 
	aa_seq2{I} = ltrs(idx) ;
	fitness_potential2(I) = sum( ltrs_fit_pot(idx) ) ;
end
idx = fitness_potential2 > 1 & fitness_potential2 < max_FP ; 
aa_seq = vertcat( aa_seq2(idx) , aa_seq) ; 
fitness_potential = vertcat( fitness_potential2(idx) , fitness_potential);  


% we sparely populate high FP values. get some more
aa_seq2 = cell( Nseqs*10 , 1);
fitness_potential2 = NaN( Nseqs*10 , 1);
for I = 1:Nseqs*10
	idx = randi(20 , genotype_length , 1); 
	aa_seq2{I} = ltrs(idx) ;
	fitness_potential2(I) = sum( ltrs_fit_pot(idx) ) ;
end
idx = fitness_potential2 > 4 & fitness_potential2 < max_FP ; 
aa_seq = vertcat( aa_seq2(idx) , aa_seq) ; 
fitness_potential = vertcat( fitness_potential2(idx) , fitness_potential); 


T = table();
T.aa_seq = aa_seq ; 
T.fitness_potential = fitness_potential ; 
T = T( T.fitness_potential > min_FP & T.fitness_potential < max_FP , : );
T = sortrows( T ,  'fitness_potential');
%% uniform distribution of fitness potentials
G = round(T.fitness_potential*10000);
T = T( diff(G)~=0 ,:) ; 

%% assign linear, sigmoid, and non-monotinc fitness-fitpot relationships

% linear
T.fitness_linear = T.fitness_potential + abs(min(T.fitness_potential)) ; 
T.fitness_linear = T.fitness_linear ./ max( T.fitness_linear );

%sigmoid1
T.fitness_logsig = logsig( T.fitness_potential - mean(abs(T.fitness_potential ))  ) ; 


% non monotonic
T.fitness_sin1 = sin( T.fitness_potential ) + abs( sin( T.fitness_potential )  ) ; 
T.fitness_sin1 = T.fitness_sin1 ./ max(T.fitness_sin1);

T.fitness_sin2 = abs(sin( T.fitness_potential ./ 2 ));

writetable( T , '~/Develop/HIS3InterspeciesEpistasis/Analysis/Lucas/Revisions/Simulations_NN_Detects_Linear_Nonmonotonic_fitness_functions.tab',...
    'FileType','text','Delimiter','\t');

fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ; 
plot( T.fitness_potential , T.fitness_sin1 , 'DisplayName' , 'sin 1' ,'LineWidth',3 )
plot( T.fitness_potential , T.fitness_sin2 , 'DisplayName' , 'sin 2' ,'LineWidth',3 )
plot( T.fitness_potential , T.fitness_linear , '-k' , 'DisplayName' , 'linear' ,'LineWidth',3  )
plot( T.fitness_potential , T.fitness_logsig , 'DisplayName' , 'sigmoid' ,'LineWidth',3  )

xlabel('Fitness potential')
ylabel('Fitness')
legend('location','best')
set(gca,'xtick',-10:10)

fh = figure('units','centimeters','position',[5 5 8 8]);
histogram( T.fitness_potential , 1e2 );
xlabel('Fitness potential')
set(gca,'xtick',-10:10)
ylabel('# of genotypes')