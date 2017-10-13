%% show distribution for a few AA seqs w/multiple synonymous variants
%  ie:  show that, for cherry picked data, our data are awesome. 

%% load synonymous variant data
T = readtable('~/Develop/HIS3InterspeciesEpistasis/Data/synonymous_variants_rescaled_data.tab','Delimiter','\t','FileType','text');

G = grpstats( T(T.SegN==7,:) , {'SegN' 'aa_seq'} ,{'min' 'max' 'mean' 'mode' 'std'} , 'DataVars','s');
G = G(G.GroupCount>8,:);
G = sortrows(G,'std_s','descend');

%%
data_fit = T.s( strcmp(T.aa_seq,'LHALAKHAGWSLIVECIGDLDIDDHHTIED') );
data_unfit = T.s( strcmp(T.aa_seq,'NNALAKHGGWSLIVECIGDLLIDDHHTAED') );
data_midfit = T.s( strcmp(T.aa_seq,'NTALAKHAGWSLIVECIGDLDIDDHHTLED') );
data_midfit2 = T.s( strcmp(T.aa_seq,'HNALAKHCGWSLIVECIGDLHIDDHHSVED') );

%%
xl = 0:0.05:1 ;
m = 1/0.47 ; 
figure;
hold on; 
histogram( m*data_fit , xl)
histogram( m*data_unfit , xl)
histogram( m*data_midfit , xl)
histogram( m*data_midfit2 , xl)
xlabel('Fitness')
ylabel('# of synonymous variants')
set(gca,'xtick',0:0.2:1)

%%  Fig #2  
%% Plot the distribution of fitness effects for some substitutions
%%  
load('/Users/lcarey/Desktop/HIS3scratch/SignEpi/SignEpi_Ronly.mat');
R = s(2).R ; 
up  = unique(R.Perm);

%%
figure; 
hold on; 
for I = 1:numel(up)
    data = R.FitImpact( strcmp(R.Perm,up{I}));
    fprintf('%d\t%s\t%0.02f\t%0.02f\n' ,I, up{I} , var(data) , var(abs(data)) );
end
%%
up  = unique(R.Perm);
figure; 
hold on; 
for I = 1:numel(up)
    data = R.FitImpact( strcmp(R.Perm,up{I}));
    fprintf('%d\t%s\t%0.02f\t%0.02f\n' ,I, up{I} , 100*mean(data>0.4) ,  100*mean(data<-0.4) );
    [f,x] = ecdf(data);
     plot(x,f,'DisplayName',up{I}); 
end
legend('location','best')

%%
fh = figure('units','centimeters','position',[5 5 7 7 ]);
hold on; 
for I = 1
    data = R.FitImpact( strcmp(R.Perm,up{I}));
    %[f,x] = ecdf(data);
    % plot(x,f,'DisplayName',up{I}); 
    histogram(data,50)
end
ylabel('Fraction of genetic backgrounds')
xlabel('Fitness impact of substitution')


print('-dpsc2','test.eps');
%legend('location','best')

%%
I = 1 ;
R2 = sortrows(R,'FitImpact');
R2 = R2(strcmp(R2.Perm,up{I}) , :);
R2 = R2( R2.Fit1<1 & R2.Fit2<1 , :);
%R2 = R2( abs(R2.FitImpact)<0.8 , :);

%%

figure; hold on ;
for I = randsample(nrows(R2),10)
    plot( [0 1 ] , [ R2.Fit1(I)  R2.Fit2(I)] ,'.-k')
end

for I = randsample(find(R2.Fit1>0.1 & R2.Fit1<0.9),5)
    plot( [0 1 ] , [ R2.Fit1(I)  R2.Fit2(I)] ,'.-k')
end

for I = randsample(find(R2.Fit2>0.1 & R2.Fit2<0.9),5)
    plot( [0 1 ] , [ R2.Fit1(I)  R2.Fit2(I)] ,'.-k')
end

xlim([-0.1 1.1])
set(gca,'xtick',0:1)
set(gca,'xticklabel',{'A' 'C'})
xlabel('Substitution')
ylim([-0.015 1.015])
ylabel('Fitness of each backgroud')
