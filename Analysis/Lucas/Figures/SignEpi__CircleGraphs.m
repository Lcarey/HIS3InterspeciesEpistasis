load('SignEpi_Pairs_04.mat');
R.ID1 = arrayfun( @(I) sprintf( '%s%d%s' , R.Perm{I}(1) , R.VarPos(I) , R.Perm{I}(2)) , 1:length(R),'UniformOutput',false)' ;
R.ID2 = arrayfun( @(I) sprintf( '%s%d%s' , R.SubPerm{I}(1) , R.SubPos(I) , R.SubPerm{I}(2)) , 1:length(R),'UniformOutput',false)' ;

allids = unique(vertcat(R.ID1,R.ID2));

%% calculate frequency of sign and reciprocal sign epistasis per position
all_reciprocal_sign = cell(12,1);
running_totals = NaN( 12 , 3);

%figname = '~/Dropbox/Pokusaeva17/Figures/SignEpiCircles2.eps';
%delete(figname);
for SegN = 1:12
    allids = unique(vertcat(R.ID1(R.SegN == SegN),R.ID2(R.SegN == SegN)));
    idx = R.SegN == SegN & R.pBon < 0.05   ; 
    ID1 = R.ID1(idx);
    ID2 = R.ID2(idx);
    logods = R.logodds(idx) ; 
    nnodes = numel(unique(vertcat(ID1,ID2)));
    nedges = sum(idx);
%     need_to_add_to_graph = allids( ~ismember( allids , vertcat(ID1,ID2)) ) ; 
%      for I = 1:numel(need_to_add_to_graph)
%          ID1 = vertcat( ID1 , need_to_add_to_graph{I});
%          ID2 = vertcat( ID2 , need_to_add_to_graph{I});
%          logods = vertcat( logods , 0.000001);
%      end
%     G = digraph( ID1 , ID2 , logods) ;
%     
%     fh = figure('units','centimeters','position',[5 5 10 10]);
%     ph = plot( G , 'layout' , 'circle' ,'ShowArrows','off'...
%         ,'LineWidth',0.25,'EdgeColor','k','NodeColor','k') ;
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     box off;
%     title( sprintf('# nodes = %d , # edges = %d' , nnodes , nedges ) )
%     text( -0.2 , 0 , num2str(SegN) ,'Color', 'r' ,'FontSize',60)
%     print('-dpsc2',figname,'-append');
%     close;
%     
    % by hand: 12 under sign; 8 under reciprocal sign
    varied_posisitions = unique( vertcat( R.VarPos( R.SegN==SegN ) , R.SubPos(R.SegN==SegN))) ; 
    FwdDir = strcat(ID1,ID2);
    RevDir = strcat(ID2,ID1);
    reciprocal_sign = intersect(FwdDir,RevDir)  ;
    sign_epi = unique(str2double(regexprep( ID1 , '[A-Z]',''))) ; 
    causing_sign_epi = unique(str2double(regexprep( ID2 , '[A-Z]',''))) ; 
    all_reciprocal_sign{SegN} = reciprocal_sign;
    
    reciprocal_sign_epi_pos = unique(str2double(regexprep(  cellfun( @(X)X(2:3) , reciprocal_sign,'UniformOutput',false) , '[A-Z]',''))) ;
    
    fprintf('S %d has %d/%d sites under sign and %d/%d under recip-sign\n' , SegN  ,  numel(causing_sign_epi) , numel(varied_posisitions) , numel(reciprocal_sign_epi_pos) , numel(varied_posisitions) );
    running_totals(SegN,:) = [ numel(causing_sign_epi) numel(reciprocal_sign_epi_pos) numel(varied_posisitions) ];
end
    
%%
fid  = fopen('RepSign.tab','w');
for I = 1:12
    for J = all_reciprocal_sign{I}'
        fprintf( fid , '%d\t%s\n' , I , char(J));
    end
end
fclose(fid);
  