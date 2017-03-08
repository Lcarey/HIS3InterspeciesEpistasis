
%% https://www.mathworks.com/help/stats/examples/fitting-data-with-generalized-linear-models.html#zmw57dd0e4667
% A set of car weights
weight = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
% The number of cars tested at each weight
tested = [48 42 31 34 31 21 23 23 21 16 17 21]';
% The number of cars failing the test at each weight
failed = [1 2 0 3 8 8 14 17 19 15 17 21]';
% The proportion of cars failing for each weight
proportion = failed ./ tested;

figure ; hold on; 
plot(weight,proportion,'s')
xlabel('Weight'); ylabel('Proportion');

[logitCoef,dev] = glmfit(weight,[failed tested],'binomial','link','logit');
logitFit = glmval(logitCoef,weight,'logit');
plot(weight,logitFit,'r-','DisplayName','logitfit');


[logitCoef,dev] = glmfit(weight,[failed tested],'binomial','link','identity');
logitFit = glmval(logitCoef,weight,'logit');
plot(weight,logitFit,'b-','DisplayName','identity');

xlabel('Weight'); ylabel('Proportion');
legend('location','nw')



%% all link functions
eta = -5:.1:5;
plot(eta,1 ./ (1 + exp(-eta)),'-', eta,normcdf(eta), '-', ...
     eta,1 - exp(-exp(eta)),'-', eta,exp(-exp(eta)),'-');
xlabel('Linear function of predictors'); ylabel('Predicted mean response');
legend('logit','probit','complementary log-log','log-log','location','east');


%%
% [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options) ;

%% 
k=10 ; x0 = 1 ; L = 1 ;
xinit = [  x0 k L ] ;
lb = [0 0 0] ;
ub = [100 100 100] ; 
logistic_function = @(x,xdata) x(3) ./ (1+exp( (-1.*x(2)).*(xdata-x(1)) ) ) ;

sim_ddG_x = linspace( 0.001 , 2 , 1e2) ;
sim_ddG_fitness = logistic_function(xinit , sim_ddG_x );
sim_ddG_fitness = sim_ddG_fitness .* random('normal',1,0.01,size(sim_ddG_x)) ;

sim_ddG_fitness_lin = linspace(0.001,1,numel(sim_ddG_x));

figure; hold on;  grid on ;
plot( sim_ddG_x , 1-sim_ddG_fitness ,'ok')
plot( sim_ddG_x , 1-sim_ddG_fitness_lin ,'ok')

[params_fit_lin,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(logistic_function,xinit,sim_ddG_x,sim_ddG_fitness_lin,lb,ub) ;
y_pred_lin = logistic_function(params_fit_lin , sim_ddG_x ) ;
[r2_lin , rmse] = rsquare(sim_ddG_fitness_lin,y_pred_lin) ;

[params_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(logistic_function,xinit,sim_ddG_x,sim_ddG_fitness,lb,ub) ;
y_pred = logistic_function(params_fit , sim_ddG_x ) ;
[r2 , rmse] = rsquare(sim_ddG_fitness,y_pred) ;

plot( sim_ddG_x , 1-y_pred ,'-r','DisplayName', sprintf('R^2=%0.02f',r2) );
plot( sim_ddG_x , 1-y_pred_lin ,'-b','DisplayName', sprintf('R^2=%0.02f',r2_lin) );

legend('location','best')
ylim([0 1])