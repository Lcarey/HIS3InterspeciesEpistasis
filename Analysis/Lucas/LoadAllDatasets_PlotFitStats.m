files = dir('March14Fits/*scaled_info.csv_.mat');
files = dir('./R_logistic_March14-5000-27/*.mat');

runtimes = NaN(0) ;
nvariants    = NaN(0);
lib = NaN(0);
filecat = NaN(0);
k = NaN(0);
test_rmse = NaN(0);
trainr2 = NaN(0);
testr2 = NaN(0);
ddGfit_corr = NaN(0,1);
testr2 = NaN(0);
test_lib_var = NaN(0,1);
pred_residual_var = NaN(0);
fitparams = cell(0);
for I = 1:numel(files)
    load( [files(I).folder filesep files(I).name]);
    metadata = regexp( files(I).name , '_' , 'split');
    
    runtimes = vertcat( runtimes, R.runtime);
    k = vertcat( k, R.k);
    testr2 = vertcat( testr2, R.test_r2);
    trainr2 = vertcat( trainr2, R.train_r2);
    
    nvariants = vertcat( nvariants,  repmat( numel(R.Train{1}) , 10 , 1 )  );
    lib = vertcat( lib , repmat(  str2double( metadata{ regexpcmp( metadata,'^S\d')  }(2:end)) , 10 , 1) );
    fitparams = vertcat( fitparams , repmat(  metadata(3) , 10 , 1) );

    filecat = vertcat( filecat , repmat( I , 10 , 1 )  ) ;
    corrs = triu(corr(cell2mat(R.fit_ddGvect)'),1) ; % correlation between all ddG fit values
    corrs = corrs(:); corrs = corrs(corrs~=0);
    corrs = mean(corrs);
    ddGfit_corr = vertcat( ddGfit_corr ,  repmat( corrs , 10 , 1 )  ) ;
   % figure; plot(R.fit_ddGvect{1},R.fit_ddGvect{2},'ok'); title(files(I).name);grid on; 
   
     % calc variance of each lib & remaining variance after model. plot residual variance
   test_lib_var = vertcat( test_lib_var , cellfun(@var , R.fit_test) ) ;
   mdls = arrayfun(@(X) GeneralizedLinearModel.fit( R.fit_test{X} , R.pred_fit_test{X}) , 1:10,'UniformOutput',false);
   pred_residual_var = vertcat( pred_residual_var , arrayfun( @(X) var(mdls{X}.Residuals.LinearPredictor) , 1:10)' ) ;

end
boxplot_order =  {   lib   } ;

figure; 
subplot(2,4,1)
boxplot( runtimes./60./60 , boxplot_order ,'notch','on')
ylabel('hours')
grid on ;

subplot(2,4,2)
boxplot( k , boxplot_order ,'notch','on')
ylabel('k')
grid on ;

subplot(2,4,3)
boxplot( testr2 , boxplot_order ,'notch','on')
ylabel('Test R^2')
grid on ;

subplot(2,4,4)
boxplot( ddGfit_corr ,boxplot_order,'notch','on')
ylabel('ddGfit corr')
grid on ;

subplot(2,4,5)
boxplot( test_lib_var ,boxplot_order,'notch','on')
ylabel('test lib var')
grid on ;

subplot(2,4,6)
boxplot( pred_residual_var ,boxplot_order,'notch','on')
ylabel('pred residual var')
grid on ;

subplot(2,4,7)
boxplot( test_lib_var-pred_residual_var ,boxplot_order,'notch','on')
ylabel('decr in  var')
grid on ;


subplot(2,4,8)
%X = grpstats( [k  testr2 ] , boxplot_order ,'median');  gscatter( X(:,2) , X(:,1) ,  unique(boxplot_order{1}) );
gscatter( trainr2 , k , boxplot_order)
xlabel('Train R^2'); set(gca,'xtick',0:0.2:1)
ylabel('K') ; set(gca,'ytick',0:10:100); 
grid on ;


% figure; hold on ; plot( k , trainr2 ,'pk','DisplayName','train');
%%
logistic_function = @(x0,k,L,sum_delta_Gs) L ./ (1+exp( (-1.*k).*(sum_delta_Gs-x0) ) ) ; % x0 , k , L



%% are there any always bad AAs?
I=4; files(I).name
load( [files(I).folder filesep files(I).name]);

AA_train_mat = R.aa_for_all_variants{1}( R.Train{1} , :);
fit_train_vect = R.fit_train{1}'; 
[ ~, fit_order] = sort(fit_train_vect,'descend');

figure ; 
imagesc( [ double(AA_train_mat(fit_order,:))./27 NaN( numel(fit_train_vect),1) fit_train_vect(fit_order) ])
colorbar;
figure ; 
imagesc( [ double(AA_train_mat)./27 NaN( numel(fit_train_vect),1) fit_train_vect ])
colorbar;
%% show fit ddG
a = R.fit_ddGvect{1};
a( round(a*10)==5 ) = NaN;
a = reshape( a , 27 , []) ; 

figure; 
minv = min(min(a));
maxv = max(max(a));
a(isnan(a)) = minv-((maxv-minv)/5);
ddd=[0.85 0.85 0.85 ;jet(10)];
colormap(ddd);
imagesc(a)
colorbar
xlabel('Position')
ylabel('Variant')
%%
MapAA2I = containers.Map(  arrayfun(@(X){X},['A':'Z' '_' ])  , 1:27 ) ; % all AAs + stop
MapI2AA = containers.Map( 1:27 ,   arrayfun(@(X){X},['A':'Z' '_' ])  ) ; % all AAs + stop

for I = 1:(size(AA_train_mat,2))
    aa_this_pos_vect = AA_train_mat(:,I);
    if(numel(unique(aa_this_pos_vect))>1)
    figure; 
    boxplot( fit_train_vect , { aa_this_pos_vect arrayfun(@(X)MapI2AA(X),aa_this_pos_vect) },'notch','on');
    title(I)
    grid on;
    ylabel('Fitness')
    end
end
%%  





