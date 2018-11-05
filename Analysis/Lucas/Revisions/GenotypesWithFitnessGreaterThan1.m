%% load data
DATADIR = '~/Develop/HIS3InterspeciesEpistasis/Data/' ; 
for I = 1:12
    s(I).T = readtable([DATADIR 'S' num2str(I) '_scaled_info_v2.csv'],'Delimiter','\t','FileType','text'); 
end
%%
figure; hold on ;
Y = NaN(0);
Y2 = NaN(0);
Y3 = NaN(0);

xl = 1:0.001:1.2 ; 
for J = xl
    M = NaN(12,1);
    M2 = NaN(12,1);
    M3 = NaN(12,1);
for I = 1:12
    data = s(I).T.s ;
    data2 =  s(I).T.s( s(I).T.size>=5) ; 
    data3 =  s(I).T.s( s(I).T.size>=10) ; 
    M(I) =  nanmean( data > J)  ;
    M2(I) =  nanmean( data2 > J)  ;
    M3(I) =  nanmean( data3 > J)  ;
end
    Y(end+1) = 100*nanmedian(M) ; 
    Y2(end+1) = 100*nanmedian(M2) ; 
    Y3(end+1) = 100*nanmedian(M3) ; 


end
plot( xl , Y , '-k','LineWidth',2)
plot( xl , Y2 , '-','LineWidth',2)
plot( xl , Y3 , '-','LineWidth',2)
ylabel('% of genotypes with fitness > X')
xlabel('Fitness threshold')
set(gca,'ytick',0:100)
set(gca,'xtick',1:0.05:3)
grid on 
legend({'all genotypes' '>=5 synonymous variants' '>=10 synonymous variants'},'location','ne')
axis tight ;

print('-dpng2','GenotypesWithFitnessGreaterThan1.png','-r300');