function s = EpistasisFindTransitionStates(DS , fit_cutoff)
%%  s = EpistasisFindTransitionStates(DS)
% 
% LBC 10/14


fit_cutoff = -2 ;
 
fast_seqs =  DS.aa_seq_short( DS.fitness >  fit_cutoff ) ;
slow_seqs =  DS.aa_seq_short( DS.fitness <= fit_cutoff ) ;

%% interesting, but not necessary
fast_seqs_hd_mat = NaN( length(fast_seqs));
for I = 1:length(fast_seqs)
    for J = I+1:length(fast_seqs)
        fast_seqs_hd_mat(I,J) = HammingDistance( fast_seqs{I},fast_seqs{J});
    end
end

figure;
imagesc(fast_seqs_hd_mat,[0,10]);colorbar;

%%
s = struct();
c=0;
for I = 1:length(fast_seqs)
    for J = I+1:length(fast_seqs)
        if HammingDistance( fast_seqs{I},fast_seqs{J})>1
            c=c+1;
      %      s(c).seq_I = fast_seqs{I};    %takes up too much memmory
      %      s(c).seq_J = fast_seqs{J};   %takes up too much memmory
            s(c).dist = fast_seqs_hd_mat(I,J);
            all_transitions = ExpandSeqAlign(fast_seqs{I},fast_seqs{J});
     %       s(c).all_transitions = all_transitions;  %takes up too much memmory
            s(c).n_transitions = length( all_transitions);
            nfast=zeros(1,length( all_transitions));nslow=nfast;ntotal=nfast;
            for PI = 1:s(c).n_transitions
                if (~strcmp(all_transitions{PI},fast_seqs{I}) & ~strcmp(all_transitions{PI},fast_seqs{J}))
                nfast(PI) = sum(ismember(fast_seqs,all_transitions{PI}));
                nslow(PI) = sum(ismember(slow_seqs,all_transitions{PI}));
                ntotal(PI) =  nslow(PI) +  nfast(PI) ;
                end
            end
            s(c).Nfast = sum( nfast  );
            s(c).Nslow = sum( nslow  );
            s(c).Nmeasured     = sum( ntotal  );
        end
    end

end


%%
%%
idx = [s.dist]>0;


figure;hold on;
[f,x]= ecdf([s(idx).Nslow] );
plot(x,f,'-b','displayname','# blocked transitions','linewidth',4);
[f,x]= ecdf([s(idx).Nfast] );
plot(x,f,'-r','displayname','# allowed transitions','linewidth',4);
[f,x]= ecdf([s(idx).Nmeasured] );
plot(x,f,'--k','displayname','# transitions measured','linewidth',1.5);

grid on;
%set(gca,'yscale','log')
legend('location','best');            
set(gca,'fontsize',15);
ylabel('fraction of AA seq pairs')
xlabel('number of transitions')
set(gca,'ytick',0:0.05:1)
set(gca,'xtick',0:10:1e4)
xlim([-2 max(xlim)])
print('-dpsc2','transitions.eps','-append');
            


figure;hold on;

[f,x]= ecdf([s(idx).Nslow]./[s(idx).N]);
plot(x,f,'-b','displayname','fraction blocked transitions','linewidth',2);
[f,x]= ecdf([s(idx).Nfast]./[s(idx).N]);
plot(x,f,'-r','displayname','fraction allowed transitions','linewidth',2);

grid on;
%set(gca,'yscale','log')
legend('location','best');            
set(gca,'fontsize',10);
ylabel('fraction of AA seq pairs')
xlabel('fraction of transitions')
set(gca,'ytick',0:0.05:1)
set(gca,'xtick',0:.05:1)
xlim([-0.01 max(xlim)])
print('-dpsc2','transitions.eps','-append');


%%
figure ; 
xl = [ 0 1:10:100 ];
c=NaN(3,length(xl));
c(1,:) = hist([s.Nmeasured],xl);
c(2,:) = hist([s.Nfast],xl);
c(3,:) = hist([s.Nslow],xl);
%%
bar(c');
set(gca,'xtick',1:length(c));
set(gca,'xticklabel',[0:10:100 ] );
legend({ 'fast transitions' ,'all found transitions' , 'slow transitions'},'location','ne')

end
