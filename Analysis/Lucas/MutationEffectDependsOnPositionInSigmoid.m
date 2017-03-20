%% The effect of mutations is context dependent based on where the starting point is in the fitness potential sigmoid

DataDir = 'Downloads/HIS3scratch/27/';
load([DataDir 'R_logistic_317-20000-27-Init05_17-Mar-201721:24:59_4315_S7_scaled_info.csv_.mat'] );

%% Pick three

figure;  hold on; 
dscatter(1-R.train_sumddG{1}' , R.fit_train{1}');
grid on ;

idx_low_ddG = find( 1-R.train_sumddG{1}' > 0.1 & 1-R.train_sumddG{1}' < 0.15 ,1,'first')
idx_mid_ddG = find( 1-R.train_sumddG{1}' > 0.26 & 1-R.train_sumddG{1}' < 0.28 ,1,'first')
idx_high_ddG = find( 1-R.train_sumddG{1}' > 0.31 & 1-R.train_sumddG{1}' < 0.4 ,1,'first')

plot( 1-R.train_sumddG{1}(idx_low_ddG) , R.fit_train{1}(idx_low_ddG),'pr');
plot( 1-R.train_sumddG{1}(idx_mid_ddG) , R.fit_train{1}(idx_mid_ddG),'pr');
plot( 1-R.train_sumddG{1}(idx_high_ddG) , R.fit_train{1}(idx_high_ddG),'pr');


low_aa = R.aa_for_all_variants{1}( R.Train{1}(idx_low_ddG) ,:) ;
mid_aa = R.aa_for_all_variants{1}( R.Train{1}(idx_mid_ddG) ,:) ;
high_aa = R.aa_for_all_variants{1}( R.Train{1}(idx_high_ddG) ,:) ;

positions_that_can_vary = find(all( vertcat(low_aa==mid_aa , low_aa==high_aa , high_aa==mid_aa ) , 1));

penaly_matrix = reshape(  R.fit_ddGvect{1}  , 27 , []);

for I = 1:size(penaly_matrix,2)
    