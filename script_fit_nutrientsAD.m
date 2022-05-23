function nuttablestat=script_fit_nutrientsAD(nuttablestat, fooddatabase0910selclean, muf,sf)


vTrunc=100;

nuttablestat.muf=muf;
nuttablestat.sf=sf;

nuttablestat.np=zeros(height(nuttablestat),1);

nuttablestat.muL=zeros(height(nuttablestat),1);
nuttablestat.sigmaL=zeros(height(nuttablestat),1);
nuttablestat.muLt=zeros(height(nuttablestat),1);
nuttablestat.sigmaLt=zeros(height(nuttablestat),1);
nuttablestat.mL=zeros(height(nuttablestat),1);
nuttablestat.sL=zeros(height(nuttablestat),1);
nuttablestat.pL=zeros(height(nuttablestat),1);
nuttablestat.ksstatL=zeros(height(nuttablestat),1);

nuttablestat.aG=zeros(height(nuttablestat),1);
nuttablestat.bG=zeros(height(nuttablestat),1);
nuttablestat.aGt=zeros(height(nuttablestat),1);
nuttablestat.bGt=zeros(height(nuttablestat),1);
nuttablestat.mG=zeros(height(nuttablestat),1);
nuttablestat.sG=zeros(height(nuttablestat),1);
nuttablestat.pG=zeros(height(nuttablestat),1);
nuttablestat.ksstatG=zeros(height(nuttablestat),1);


nuttablestat.aW=zeros(height(nuttablestat),1);
nuttablestat.bW=zeros(height(nuttablestat),1);
nuttablestat.aWt=zeros(height(nuttablestat),1);
nuttablestat.bWt=zeros(height(nuttablestat),1);
nuttablestat.mW=zeros(height(nuttablestat),1);
nuttablestat.sW=zeros(height(nuttablestat),1);
nuttablestat.pW=zeros(height(nuttablestat),1);
nuttablestat.ksstatW=zeros(height(nuttablestat),1);



nuttablestat.muNt=zeros(height(nuttablestat),1);
nuttablestat.sigmaNt=zeros(height(nuttablestat),1);
nuttablestat.mN=zeros(height(nuttablestat),1);
nuttablestat.sN=zeros(height(nuttablestat),1);
nuttablestat.pN=zeros(height(nuttablestat),1);
nuttablestat.ksstatN=zeros(height(nuttablestat),1);


for ind=1:size(fooddatabase0910selclean,1)
    disp('****************************************************************')
    disp(ind)
    v=fooddatabase0910selclean(ind,:);
    disp(sum(v>0))
    disp([nuttablestat.muf(ind), nuttablestat.sf(ind)])
    v=v(v>0);
    nuttablestat.np(ind)=length(v);
    if or(or(or(nuttablestat.muf(ind)==0, nuttablestat.sf(ind)==0), isnan(nuttablestat.muf(ind))), isnan(nuttablestat.sf(ind)))
        disp(ind)
        disp("No content")
        nuttablestat.muLt(ind)=nan;
        nuttablestat.sigmaLt(ind)=nan; 
        nuttablestat.mL(ind)=nan;
        nuttablestat.sL(ind)=nan;
        nuttablestat.pL(ind)=nan;
        nuttablestat.ksstatL(ind)=nan; 
        
        nuttablestat.aGt(ind)=nan;
        nuttablestat.bGt(ind)=nan; 
        nuttablestat.mG(ind)=nan;
        nuttablestat.sG(ind)=nan;
        nuttablestat.pG(ind)=nan;
        nuttablestat.ksstatG(ind)=nan;  
        

        nuttablestat.aWt(ind)=nan;
        nuttablestat.bWt(ind)=nan; 
        nuttablestat.mW(ind)=nan;
        nuttablestat.sW(ind)=nan;
        nuttablestat.pW(ind)=nan;
        nuttablestat.ksstatW(ind)=nan; 
        
        nuttablestat.muNt(ind)=nan;
        nuttablestat.sigmaNt(ind)=nan; 
        nuttablestat.mN(ind)=nan;
        nuttablestat.sN(ind)=nan;
        nuttablestat.pN(ind)=nan;
        nuttablestat.ksstatN(ind)=nan;     
        
        
    else   
        
    % lognormal
    pdf_truncln = @(v,mu,sigma) lognpdf(v,mu,sigma) ./ logncdf(vTrunc,mu,sigma);
    [lognp,lognci] = lognfit(v);
    start = [lognp(1),lognp(2)];
    nuttablestat.muL(ind)=lognp(1);
    nuttablestat.sigmaL(ind)=lognp(2);
    try
    [paramEstsL,paramCIsL] = mle(v, 'pdf',pdf_truncln, 'start',start, 'lower',[-Inf 0]);
    nuttablestat.muLt(ind)=paramEstsL(1);
    nuttablestat.sigmaLt(ind)=paramEstsL(2);  
    pd = makedist('Lognormal', paramEstsL(1), paramEstsL(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    nuttablestat.mL(ind)=mean(r);
    nuttablestat.sL(ind)=std(r);
    [~,pL,ksstatL,~]= kstest(v,'CDF',t);
    nuttablestat.pL(ind)=pL;
    nuttablestat.ksstatL(ind)=ksstatL;
    catch
        disp(ind)
        disp("Truncated Lognormal")
        disp("No convergence")
        disp(nuttablestat.np(ind))
        
        nuttablestat.muLt(ind)=nan;
        nuttablestat.sigmaLt(ind)=nan; 
        nuttablestat.mL(ind)=nan;
        nuttablestat.sL(ind)=nan;
        nuttablestat.pL(ind)=nan;
        nuttablestat.ksstatL(ind)=nan;        
    end   
        
 

    
    
    % gamma
    pdf_truncgamma = @(v,a,b) gampdf(v,a,b) ./ gamcdf(vTrunc,a,b);
    [gammap,gammaci] = gamfit(v);
    start = [gammap(1),gammap(2)];
    nuttablestat.aG(ind)=gammap(1);
    nuttablestat.bG(ind)=gammap(2);
    try
    [paramEstsG,paramCIsG] = mle(v, 'pdf',pdf_truncgamma, 'start',start, 'lower',[0 0]);
    nuttablestat.aGt(ind)=paramEstsG(1);
    nuttablestat.bGt(ind)=paramEstsG(2);  
    pd = makedist('Gamma', paramEstsG(1), paramEstsG(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    nuttablestat.mG(ind)=mean(r);
    nuttablestat.sG(ind)=std(r);
    [~,pG,ksstatG,~]= kstest(v,'CDF',t);
    nuttablestat.pG(ind)=pG;
    nuttablestat.ksstatG(ind)=ksstatG;
    catch
        disp(ind)
        disp("Truncated Gamma")
        disp("No convergence")
        disp(nuttablestat.np(ind))
        
        nuttablestat.aGt(ind)=nan;
        nuttablestat.bGt(ind)=nan; 
        nuttablestat.mG(ind)=nan;
        nuttablestat.sG(ind)=nan;
        nuttablestat.pG(ind)=nan;
        nuttablestat.ksstatG(ind)=nan;     
    end
    
 
    % weibull
    pdf_truncwbl = @(v,a,b) wblpdf(v,a,b) ./ wblcdf(vTrunc,a,b);
    [wblp,wblci] = wblfit(v);
    start = [wblp(1),wblp(2)];
    nuttablestat.aW(ind)=wblp(1);
    nuttablestat.bW(ind)=wblp(2);
    try
    [paramEstsW,paramCIsW] = mle(v, 'pdf',pdf_truncwbl, 'start',start, 'lower',[0 0]);
    nuttablestat.aWt(ind)=paramEstsW(1);
    nuttablestat.bWt(ind)=paramEstsW(2);  
    pd = makedist('Weibull', paramEstsW(1), paramEstsW(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    nuttablestat.mW(ind)=mean(r);
    nuttablestat.sW(ind)=std(r);
    [~,pW,ksstatW,~]= kstest(v,'CDF',t);
    nuttablestat.pW(ind)=pW;
    nuttablestat.ksstatW(ind)=ksstatW;
    catch
        disp(ind)
        disp("Truncated Weibull")
        disp("No convergence")
        disp(nuttablestat.np(ind))
        
        nuttablestat.aWt(ind)=nan;
        nuttablestat.bWt(ind)=nan; 
        nuttablestat.mW(ind)=nan;
        nuttablestat.sW(ind)=nan;
        nuttablestat.pW(ind)=nan;
        nuttablestat.ksstatW(ind)=nan;     
    end   
    

    % gaussian
    pdf_truncnorm = @(v,mu,sigma) normpdf(v,mu,sigma) ./ (normcdf(vTrunc,mu,sigma)-normcdf(0,mu,sigma));
    start = [mean(v),std(v)];
    try
    [paramEstsN,paramCIsN] = mle(v, 'pdf',pdf_truncnorm, 'start',start, 'lower',[0 0]);
    nuttablestat.muNt(ind)=paramEstsN(1);
    nuttablestat.sigmaNt(ind)=paramEstsN(2);   
    pd = makedist('Normal', paramEstsN(1), paramEstsN(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    nuttablestat.mN(ind)=mean(r);
    nuttablestat.sN(ind)=std(r);    
    [~,pN,ksstatN,~]= kstest(v,'CDF',t);
    nuttablestat.pN(ind)=pN;
    nuttablestat.ksstatN(ind)=ksstatN;    
    catch
        disp(ind)
        disp("Truncated Gaussian")
        disp("No convergence")
        disp("Fit Half Gaussian")
        disp(nuttablestat.np(ind))
        
        
        
        pd = fitdist(v','HalfNormal');
        nuttablestat.muNt(ind)=pd.mu;
        nuttablestat.sigmaNt(ind)=pd.sigma;
        r = random(pd,length(v),1);
        nuttablestat.mN(ind)=mean(r);
        nuttablestat.sN(ind)=std(r);    
        [~,pN,ksstatN,~]= kstest(v,'CDF',pd);
        nuttablestat.pN(ind)=pN;
        nuttablestat.ksstatN(ind)=ksstatN; 
        
    end
    
    
    end
end