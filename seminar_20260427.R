

################## generate artificial data of new infections ##################
rm(list = ls()) 
library(deSolve);
derivs=function(t,state,pars){
      with(as.list(c(state,pars)),{
 
      if(t<=50){           
        c_t = 10;
      }else if(t<=200){
        c_t = 10-0.05*(t-100)
         }else{
        c_t = 5 + 0.03*(t-200);
       }

            dS=-beta*c_t*S*I/(S+E+I+R)
            dE=beta*c_t*S*I/(S+E+I+R)-sigma*E
            dI=sigma*E-gamma*I
            dR=gamma*I
           
            dnew = beta*c_t*S*I/(S+E+I+R);
            dcon = gamma*I

            return(list(c(dS,dE,dI,dR,dnew,dcon)))
          }
          )
    }  

beta=0.055; sigma=1/4; gamma=1/4;
R0=10*beta/gamma
print(paste0('R0=',R0))
y0=c(S=20000000,E=50,I=10,R=0,dnew=50,dcon=0);
parameters=c(beta,sigma,gamma);
sample=ode(y0,times=1:110,derivs,parameters);

tmp = sample[,6];
new = 1:360; new[1] = tmp[1];
new[2:360] = tmp[2:360]-tmp[2:360-1];

onset=sample[,3]*sigma;
data=rep(0,length(onset));
for(i in 1:length(onset)){data[i]=rpois(1,onset[i])}
################## generate artificial data of new infections ##################


################## estimate R0 using exponential growth method ##################
y=data[1:50];x=1:50;
r=coefficients(glm(y~x,family=poisson()) )[2];

E=8; V=32;
b=V/E;a=E/b;
f<-function(x){exp(-r*x)*dgamma(x,shape=a,scale=b)}
R0.m1=1/integrate(f,0,50)$value; 
print(paste0('estimated R0 using exponential growth method: ',R0.m1))
################## estimate R0 using exponential growth method ##################


################## estimate R0 using renewal equation ##################
g=1:360; g[1]=pgamma(1.5,shape=a,scale=b);
for(i in 2:length(g)){g[i]=pgamma(i+0.5,shape=a,scale=b)-pgamma(i-0.5,shape=a,scale=b)}

L<-function(R0){
    likelihood=rep(0,50);lambda=rep(0,50);
    for(i in 1:50){
      lambda[i]=R0*sum(data[1:(i-1)]*g[(i-1):1])}

    for(i in 1:50){
      likelihood[i]=data[i]*log(lambda[i])-lambda[i]}

return(sum(likelihood));
}
R0.m2=optimize(L,lower=1,upper=10,maximum=T)$maximum;
print(paste0('estimated R0 using renewal equation method: ',R0.m2))
################## estimate R0 using renewal equation ##################


################## estimate Rt using renewal equation ##################
Rt.m1=rep(0,length(data));
for(i in 2:length(data)){Rt.m1[i]=data[i]/sum(data[1:(i-1)]*g[(i-1):1])}
################## estimate Rt using renewal equation ##################


################## estimate Rt using Epiestim method ##################
Rt.m2=rep(0,length(data));
tau=7;
b=25/5;a=5/b;
for(i in tau:length(data)){
  shape_p=a+sum(data[(i-tau+1):i]);
  
  s=rep(0,tau);
  for(j in (i-tau+1):i){s[j-(i-tau)]=sum(data[1:(j-1)]*g[(j-1):1])}
  scale_p=1/(1/b+sum(s));

  Rt.m2[i]=shape_p*scale_p;
}
Rt.m2[1:(tau-1)]=Rt.m2[tau];

plot(Rt.m1,xaxs="i",yaxs='i',col=1,pch='',ylim=c(0,5));
lines(Rt.m1,col=1,lwd=2);
lines(Rt.m2,col=2,lwd=2);
abline(h=1,lty=2);
legend('topright',cex=2,inset=c(0.12,0.015),legend=c("method 1","method 2"),col=c(1,2),pch=c('',''),lwd=2,lty=1,bg='white');
points(data/max(data)*4.8);
################## estimate Rt using Epiestim method ##################


################## estimate Rt using Epiestim package ##################
install.packages('EpiEstim');
library(EpiEstim);
Rt=estimate_R(incid =data,method = "parametric_si",config = make_config(list(mean_si =E,std_si =sqrt(V)))); 
plot(Rt);

te=50:length(data);ts=te-6;ts[1]=2;
Rt=estimate_R(incid =data,method = "parametric_si",config = make_config(list(mean_si =E,std_si =sqrt(V),t_start=ts,t_end=te))); 
plot(Rt);
################## estimate Rt using Epiestim package ##################


