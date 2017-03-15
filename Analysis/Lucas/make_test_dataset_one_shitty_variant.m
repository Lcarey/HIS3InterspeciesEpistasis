%% generate a test dataset in which one AA variant is always bad

%% load data
figure; hold on;
T = readtable('~/Develop/HIS3InterspeciesEpistasis/Data/S7_scaled_info.csv'); ecdf(T.s);
%T = T(1:5000,:);

T = T( T.middle & T.nogap & ~T.stop & ~T.nonsense , :) ; ecdf(T.s);
T = T( logical(T.lib) , :); ecdf(T.s);
T = T( logical(T.nat_lib) , :); ecdf(T.s);

legend({'1' '2' '3' '4'})

%% are there any always bad AAs?
for I = 1:(mode(T.len))
    aa_vars =  cellfun( @(X)X(I) , T.aa_seq);
    if(numel(unique(aa_vars))>1)
    figure; 
    boxplot( T.s ,aa_vars);
    title(I)
    grid on;
    ylabel('Fitness')
    end
end
%%  make X(1)==Y always dead
%     every thing else always perfect
%     model should say that only one position matters

T.s = ones( height(T),1);
first_aa =  cellfun( @(X)X(2) , T.aa_seq);
T.s( first_aa=='T') = 0 ; 

writetable( T , 'S0_scaled_info.csv', 'FileType','text','Delimiter','\t')

%% did the results find that that one variant is the worse? 
T = readtable('S0_scaled_info.csv');
%%
ddGVECT = vertcat( R.fit_ddGvect{1}  , R.fit_ddGvect{2} , R.fit_ddGvect{3} , R.fit_ddGvect{4});

n_positions = T.len(1)  ;

ddGmat = reshape( mean(ddGVECT) , [], n_positions ) ;
figure; 
imagesc(ddGmat);
xlabel('Position')
ylabel('AA variant')
colorbar;
set(gca,'xtick',0:2:100)
grid on;
