nS = 15;
nR = 10000;
di = randn(1,nR,nS)+5;
C = 0.5* rand(1,nR,nS)+0.5*randn(1,nR,nS) + di;
randPerc = randn(1,nR,nS);
A = randPerc.*C;
B = (1-randPerc).*C;


% function [Ty,diff,se,CI,p,tcrit,df]=limo_yuend_ttest(a,b,percent,alpha)
% function [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alphav, nullvalue)
t2 = [];t1=[];tDi=[];
t2 = limo_yuend_ttest(A,B);
t1 = limo_trimci(A-B);
tDi =limo_trimci(di);
figure,hist([t1;t2]')