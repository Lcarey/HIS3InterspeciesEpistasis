function EpistasisLBC__SynVariantsPerSegment()
%% EpistasisLBC__SynVariantsPerSegment
% for each segment, calc various reproducibility stats
%
% LBC

DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
fns = dir( [DataDir filesep 'S*_sum_*fit.csv']);


for I = 1:12
    metadata = regexp( fns(I).name , '_' ,'split')
    fns(I).segment = str2double(metadata{1}(2:end));
    T = readtable( [ fns(I).folder filesep fns(I).name ] );
    G = grpstats( T , 'aa_seq' , {'median' 'mean' 'std' 'CoeffVar'} ,'DataVars','s');
    G.Segment =  repmat( fns(I).segment , height(G) , 1);
    G = G( G.GroupCount>=3 ,:);
    if I==1
        R=G;
    else
        R = vertcat(R,G);
    end
end

%%
RRR = grpstats( R.std_s , R.Segment ,'median') ;
[~,o] = sort(RRR,'descend');

clrs = parula(13);
figure; hold on ; grid on ;
for I = o'
    X = R.std_s(R.Segment==I);
    [f,x]=ecdf(X);
    plot(x,f,'-','Color',clrs(I,:),'DisplayName', num2str(I) )
end
xlim([0 2])
legend('location','se');

%%
figure ;
boxplot( R.std_s , R.Segment ,'Symbol','','Notch','on');
ylim([ 0 1])
grid on; 
ylabel('Std dev of syn variants')
xlabel('Segment')


%% 
DataDir = '~/Develop/HIS3InterspeciesEpistasis/Data/';
fns = dir( [DataDir filesep 'S*_sum_*fit.csv']);

figure; 
for I = 1:12
    subplot(2,6,I)
    metadata = regexp( fns(I).name , '_' ,'split') ;
    fns(I).segment = str2double(metadata{1}(2:end));
    T = readtable( [ fns(I).folder filesep fns(I).name ] );
    T = T( ~isnan(T.s) , :);
    G = grpstats( T , 'aa_seq' , {'median' 'mean' 'std' 'CoeffVar'} ,'DataVars','s');
    T.cons_s = arrayfun( @(X)G.mean_s( ismember(G.aa_seq,T.aa_seq{X})) , 1:height(T))' ;
    
end