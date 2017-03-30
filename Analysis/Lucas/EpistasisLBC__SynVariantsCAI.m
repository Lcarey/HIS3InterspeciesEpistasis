function s = EpistasisLBC__SynVariantsCAI(SegI , remove_outliers_method,remove_outliers_value )
%% Do synonymous variants that differ by codon bias have different expression?
%% EpistasisLBC__SynVariantsCAI.m
%% load data
%%

function_start_time = char(datetime) ;
if ~exist('remove_outliers_method','var')
	remove_outliers_method='std';
end
if ~exist('remove_outliers_value','var')
	remove_outliers_value=2;
end


DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/Synonymous/';
fns = dir( [DataDir filesep 'S*_sum_*fit.csv'])
fns(SegI)
%for SegI = 1:12
metadata = regexp( fns(SegI).name , '_' ,'split');
fns(SegI).segment = str2double(metadata{1}(2:end));
source_file_name = [ DataDir filesep fns(SegI).name ] ;


output_mat_file_name  = sprintf('SynVariantsCAI_%02d_%s%d_%s.mat' , fns(SegI).segment , remove_outliers_method , remove_outliers_value ,  function_start_time)

T = readtable( source_file_name );
T = T( cellfun( @length , T.seq) ==  mode(cellfun( @length , T.seq)) , :); % only mode seq length
T = T( ~isnan(T.s) , :);
% Calc CAI / nTE , etc
tic;
T.nTE  = cellfun(@(X)calc_tAI(X,'nTE'),T.seq) ;
T.tAI  = cellfun(@(X)calc_tAI(X,'tAI'),T.seq) ;
fprintf('Calc CAI took %0.01f sec\n' , toc);
% Correlation between CAI & fitness in different ranges of mean AA fitness?

tic;
G = grpstats(T ,'aa_seq' ,'median','DataVars','s');
G = G( G.GroupCount  > 1 , :);
D = T( ismember( T.aa_seq , G.aa_seq) ,:);
C = NaN( numel(G.aa_seq) , 8);

nTE_all = cell( height(G) , 1);
tAI_all = cell( height(G) , 1);
fit_all = cell( height(G) , 1);

keep_rows_idx = true(height(G),1);
for I = 1:height(G)
    idx = find(strcmp( D.aa_seq , G.aa_seq{I}));
    [s,o] = sort(D.s(idx)) ;
    nTE =  D.nTE(idx(o));
    tAI = D.tAI(idx(o));

	if strcmp(remove_outliers_method,'std')
		KeepIDX = s > ( median(s)-remove_outliers_value*std(s)) & s < ( median(s)+remove_outliers_value*std(s)) ;
		s = s(KeepIDX);
		nTE = nTE(KeepIDX);
		tAI = tAI(KeepIDX);
	elseif strcmp(remove_outliers_method,'abs')
		KeepIDX = s > ( median(s)-remove_outliers_value) & s < ( median(s)+remove_outliers_value) ;
		s = s(KeepIDX);
		nTE = nTE(KeepIDX);
		tAI = tAI(KeepIDX);
	end
    nTE_all{I} = nTE;
    tAI_all{I} = tAI;
    fit_all{I} = s;

	if numel(tAI)<2
		keep_rows_idx(I)=false;
	else

    if numel(tAI)>=6
        C(I,3) = log2( mean(nTE(end-2:end)) /  mean(nTE(1:3)) );
        C(I,6) = log2( mean(tAI(end-2:end)) /  mean(tAI(1:3)) );
    end
    C(I,7) = log2( mean(nTE(end)) /  mean(nTE(1)) );
    C(I,8) = log2( mean(tAI(end)) /  mean(tAI(1)) );
    
    [ C(I,1) , C(I,2) ] = corr(s,nTE) ;
    [ C(I,4) , C(I,5) ] = corr(s,tAI) ;
	end
    
end
G.nTE_C = C(:,1);
G.nTE_P = C(:,2);
G.tAI_C = C(:,4);
G.tAI_P = C(:,5);
G.nTE_l2_3 = C(:,3);
G.tAI_l2_3 = C(:,6);
G.nTE_l2_1 = C(:,7);
G.tAI_l2_1 = C(:,8);

G.tAI_all = tAI_all ;
G.fit_all = fit_all ;
G.nTE_all = nTE_all ;
G.segment = repmat( uint8(fns(SegI).segment) , height(G) , 1);
G = G( keep_rows_idx , :);
fprintf('Calc corr took %0.01f sec\n' , toc);

fns(SegI).T = T ;
fns(SegI).G = G ;

SegI
s = fns(SegI) ;
s.source_file_name = source_file_name ;
s.remove_outliers_method = remove_outliers_method ; 
s.remove_outliers_value  = remove_outliers_value ;
save( output_mat_file_name , 's' ) ;

T(1:10,:)
G(1:10,:)

