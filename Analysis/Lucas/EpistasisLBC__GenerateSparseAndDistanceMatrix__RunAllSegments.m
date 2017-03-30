%% EpistasisLBC__GenerateSparseAndDistanceMatrix__RunAllSegments

%% how does each matrix correlate w/fitness? 
R  = EpistasisLBC__GenerateSparseAndDistanceMatrix(1);
R.N = repmat(1,size(R,1),1)


for I = 2:12
    Q  = EpistasisLBC__GenerateSparseAndDistanceMatrix(I);
    Q.N  = repmat(I,size(Q,1),1);
    R = vertcat(Q,R);
    I
end
save('scoring_matricies.mat','R')

%%  normalize each segment to find which matricies do the best
Gs = grpstats(R,'N','mean','DataVars','corrs') ;
Gm = grpstats(R,'matricies','mean','DataVars','corrs') ; 

R.corr_avg_this_segment = NaN(length(R),1);
R.corr_avg_this_matrix  = NaN(length(R),1);

for I = 1:length(R)
    this_matrix = R.matricies{I};
    this_segment = R.N(I);
    
    R.corr_avg_this_matrix(I)  = Gm.mean_corrs( strcmp( Gm.matricies , this_matrix)); 
    R.corr_avg_this_segment(I) = Gs.mean_corrs( Gs.N == this_segment) ;
end

R.corr_norm_by_segment = R.corrs ./ R.corr_avg_this_segment  ;
R.corr_norm_by_matrix  = R.corrs ./ R.corr_avg_this_matrix  ;

%%
figure; 
plot(  R.corrs( strcmp(R.matricies,'DAYHOFF'))  ,  R.corrs( strcmp(R.matricies,'GONNET'))  ,'ok');
xlabel('DAYHOFF')
ylabel('GONNET')
xlim([-.1 1]),ylim(xlim),line(xlim,xlim)

%%
figure; 
boxplot(R.corr_norm_by_segment , R.matricies )
grid on ;

%%
figure; 
boxplot(R.corr_norm_by_matrix , R.N )
grid on ;