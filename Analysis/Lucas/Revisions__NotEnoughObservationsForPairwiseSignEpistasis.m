% LBC September 1, 2017
% d) To handle epistasis, the authors suddenly switch mathematical frameworks, employing scores based on fractal dimensionality of interacting variant pairs. This discussion suffers from a number of problems. First, it does not appear to be discussed mathematically or statistically in the methods section (though perhaps in the missing lines between 452 and 490). Second, what is its relationship to the original deep learning framework?
% An integrated method of handling epistatic interactions would be to add cross terms to the net in the first layer. 
% Third, the claims that the non-epistatic model is very good and that 40% of sites participate in epistatic interactions are in tension, and not very well resolved by the authors.

% There are too few data points to be able to do this. In short, there are many combinations and to do a neural-net you need at least XXX observations per each variable in the matrix, and in our case we have YYY. This is the reason why we have switched the frameworks, which we now explain the manuscript.

load('/Users/lcarey/Desktop/HIS3scratch/SignEpi/SignEpi_Pairs_04.mat')
%%
sigeffect = R.pBon < 0.05 ;
% sigeffect = true(size(sigeffect)); % to see all pairs
minor_topleft = cellfun( @(X) X(1,1) , R.X ); 
minor_bottomright = cellfun( @(X) X(2,2) , R.X ); 
major_topright = cellfun( @(X) X(1,2) , R.X ); 
major_bottomleft = cellfun( @(X) X(2,1) , R.X ); 

fh = figure('units','centimeters','position',[5 5 6 6 ]);
%histogram( vertcat( minor_topleft , minor_bottomright ) , 1:10:1e3 )
histogram( vertcat( minor_topleft(sigeffect) , minor_bottomright(sigeffect) ) , 1:25:1e3 )
ylabel('# of pairs of substitutions')
xlabel('# of genetic backgrounds')
%histogram( vertcat( minor_topleft , minor_bottomright ) , 1:10:1e3 )
[f,x] = ecdf( vertcat( minor_topleft(sigeffect) , minor_bottomright(sigeffect) ));
fh = figure('units','centimeters','position',[5 5 6 6 ]);
plot(x,f*100,'-k','LineWidth',3)
xlim([0 400])
ylim([0 100])
ylabel('% of pairs of substitutions')
xlabel('# of genetic backgrounds')
set(gca,'xtick',0:100:1e4)
grid on ;
set(gca,'ytick',0:25:100)