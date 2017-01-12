
proc import datafile="C:/Users/dcries/test.csv" out=mydata dbms=csv replace;
    getnames=yes;
run;

proc print data=mydata; run;

proc nlmixed data=mydata;

parms   a1=-1 a2=-1 a3=-1 b1=0 b2=0 b3=0 taue=1;

bounds va vb taue>=0;

teta = a1*x1+a2*x2+a3*x3 ;
expteta = exp(-teta);
p=1/(1+expteta);

if y=0 then loglik = log(1-p);

if y>0 then do;
	mu = exp(b1*x1+b2*x2+b3*x3);
	logden = -log(y*(2*3.14159)**0.5) - 0.5*log(taue) - ((log(y)-mu)**2)/(2*taue);
	/*logden = -log((2*3.14159)**0.5) - 0.5*log(taue) - ((log(y)-mu)**2)/(2*taue);*/
	logLik = log(p) + logden;
end;

model y ~ general(loglik);

run;


proc import datafile="C:/Users/dcries/github/bouts/data/FinalPAMSDataSetNew.csv" out=pams dbms=csv replace;
    getnames=yes;
run;



proc sort data=pams;
	by id;
	run;

data pams;
	set pams;
	if ModPAR=0 then YN=0;
	else if ModPAR>0 then YN=1;
run;

proc print data=pams(obs=10);run;
/*model without any random effects*/
	/*run first to get starting values*/
proc nlmixed data=pams;

parms  alpha0=0 alpha1=0 alpha2=0 alpha3=0
		beta0=0 beta1=0 beta2=0 beta3=0 taue=1

		;


bounds taue>=0;

teta = alpha0+alpha1*Age+alpha2*BMI+alpha3*Gender;
expteta = exp(-teta);
p=1/(1+expteta);

if ModPAR=0 then loglik = log(1-p);

if ModPAR>0 then do;
	mu = exp(beta0+beta1*Age+beta2*BMI+beta3*Gender);
	logden = -log(ModPAR*(2*3.14159)**0.5) - 0.5*log(taue) - ((log(ModPAR)-mu)**2)/(2*taue);
	/*logden = -log((2*3.14159)**0.5) - 0.5*log(taue) - ((log(y)-mu)**2)/(2*taue);*/
	logLik = log(p) + logden;
end;

model ModPAR ~ general(loglik);
run;

/*model with uncorrelated random effects*/
proc nlmixed data=pams;

parms  a=0 b=0 alpha0=2.28 alpha1=0 alpha2=-0.02
		beta0=1.61 beta1=0 beta2=0 taue=1.26

		var1=1 var2=1 ;


bounds var1 var2 taue>=0;

teta = a+alpha0+alpha1*Age+alpha2*BMI;
expteta = exp(-teta);
p=1/(1+expteta);

if ModPAR=0 then loglik = log(1-p);

if ModPAR>0 then do;
	mu = exp(b+beta0+beta1*Age+beta2*BMI);
	logden = -log(ModPAR*sqrt(2*3.14159)) - 0.5*log(taue) - ((log(ModPAR)-mu)**2)/(2*taue);
	/*logden = -log((2*3.14159)**0.5) - 0.5*log(taue) - ((log(y)-mu)**2)/(2*taue);*/
	logLik = log(p) + logden;
end;

model ModPAR ~ general(loglik);
random a b ~ normal([0,0],[var1,0,var2]) subject=id;
run;

proc nlmixed data=pams;

parms  a=0 b=0 alpha0=0 alpha1=0 alpha2=0
		beta0=0 beta1=0 beta2=0 taue=1

		var1=.1 var2=.1 cov12 = 0.01;


bounds var1 var2 taue>=0;

teta = a+alpha0+alpha1*Age+alpha2*BMI;
expteta = exp(-teta);
p=1/(1+expteta);

*mu = exp(b+beta0+beta1*Age+beta2*BMI);
mu = (b+beta0+beta1*Age+beta2*BMI);
llb=log((1-p)**(1-YN)) + log(p**(YN));
/*if ModPAR>0 then do; ll1=log(1/(sqrt(2*3.14159*taue)*ModPAR));
                  ll2=(-(log(ModPAR)-mu)**2)/(2*taue);
                  ll_ln=ll1+ll2; end;*/
if ModPAR>0 then do; ll1=log(1/(sqrt(2*3.14159*taue)));
                  ll2=(-(log(ModPAR)-mu)**2)/(2*taue);
                  ll_ln=ll1+ll2; end;


if ModPAR=0 then loglik=llb;
else if ModPAR>0 then loglik=llb+ll_ln;

/*if ModPAR=0 then loglik = log(1-p);

if ModPAR>0 then do;
	mu = exp(b+beta0+beta1*Age+beta2*BMI);
	logden = -log(ModPAR*sqrt(2*3.14159*taue)) - ((log(ModPAR)-mu)**2)/(2*taue);
	logLik = log(p) + logden;
end; */

model ModPAR ~ general(loglik);
random a b ~ normal([0,0],[var1,cov12,var2]) subject=id;
run;
