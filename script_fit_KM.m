function tGsel=script_fit_KM(tGsel, db)

tGsel.muL=zeros(height(tGsel),1);
tGsel.sigmaL=zeros(height(tGsel),1);
tGsel.mL=zeros(height(tGsel),1);
tGsel.sL=zeros(height(tGsel),1);
tGsel.pL=zeros(height(tGsel),1);
tGsel.ksstatL=zeros(height(tGsel),1);

tGsel.aG=zeros(height(tGsel),1);
tGsel.bG=zeros(height(tGsel),1);
tGsel.mG=zeros(height(tGsel),1);
tGsel.sG=zeros(height(tGsel),1);
tGsel.pG=zeros(height(tGsel),1);
tGsel.ksstatG=zeros(height(tGsel),1);


tGsel.aW=zeros(height(tGsel),1);
tGsel.bW=zeros(height(tGsel),1);
tGsel.mW=zeros(height(tGsel),1);
tGsel.sW=zeros(height(tGsel),1);
tGsel.pW=zeros(height(tGsel),1);
tGsel.ksstatW=zeros(height(tGsel),1);



tGsel.muN=zeros(height(tGsel),1);
tGsel.sigmaN=zeros(height(tGsel),1);
tGsel.mN=zeros(height(tGsel),1);
tGsel.sN=zeros(height(tGsel),1);
tGsel.pN=zeros(height(tGsel),1);
tGsel.ksstatN=zeros(height(tGsel),1);



tGsel.muE=zeros(height(tGsel),1);
tGsel.mE=zeros(height(tGsel),1);
tGsel.sE=zeros(height(tGsel),1);
tGsel.pE=zeros(height(tGsel),1);
tGsel.ksstatE=zeros(height(tGsel),1);

tGsel.aU=zeros(height(tGsel),1);
tGsel.bU=zeros(height(tGsel),1);
tGsel.mU=zeros(height(tGsel),1);
tGsel.sU=zeros(height(tGsel),1);
tGsel.pU=zeros(height(tGsel),1);
tGsel.ksstatU=zeros(height(tGsel),1);


for ind=1:height(tGsel)
    eclabel=tGsel.EC(ind);
    sublabel=tGsel.substrate(ind);
    
    filter=and(db.EC==eclabel,db.substrate==sublabel);
    tsel=db(filter,:);
    v=tsel.KM;
    
     % lognormal
    [lognp,lognci] = lognfit(v);
    tGsel.muL(ind)=lognp(1);
    tGsel.sigmaL(ind)=lognp(2);
    pd = makedist('Lognormal', lognp(1), lognp(2));
    r = random(pd,length(v),1);
    tGsel.mL(ind)=mean(r);
    tGsel.sL(ind)=std(r);
    [~,pL,ksstatL,~]= kstest(v,'CDF',pd);
    tGsel.pL(ind)=pL;
    tGsel.ksstatL(ind)=ksstatL;

        
    
    % gamma
    [gammap,gammaci] = gamfit(v);
    tGsel.aG(ind)=gammap(1);
    tGsel.bG(ind)=gammap(2);  
    pd = makedist('Gamma', gammap(1), gammap(2));
    r = random(pd,length(v),1);
    tGsel.mG(ind)=mean(r);
    tGsel.sG(ind)=std(r);
    [~,pG,ksstatG,~]= kstest(v,'CDF',pd);
    tGsel.pG(ind)=pG;
    tGsel.ksstatG(ind)=ksstatG;

    
 
    % weibull
    [wblp,wblci] = wblfit(v);
    tGsel.aW(ind)=wblp(1);
    tGsel.bW(ind)=wblp(2);
    pd = makedist('Weibull', wblp(1), wblp(2));
    r = random(pd,length(v),1);
    tGsel.mW(ind)=mean(r);
    tGsel.sW(ind)=std(r);
    [~,pW,ksstatW,~]= kstest(v,'CDF',pd);
    tGsel.pW(ind)=pW;
    tGsel.ksstatW(ind)=ksstatW;
    
    
    %gaussian   
    pdf_truncnorm = @(v,mu,sigma) normpdf(v,mu,sigma) ./ (normcdf(inf,mu,sigma)-normcdf(0,mu,sigma));
    start = [mean(v),std(v)];
    try
    [paramEstsN,paramCIsN] = mle(v, 'pdf',pdf_truncnorm, 'start',start, 'lower',[0 0]);
    tGsel.muN(ind)=paramEstsN(1);
    tGsel.sigmaN(ind)=paramEstsN(2);   
    pd = makedist('Normal', paramEstsN(1), paramEstsN(2));
    t = truncate(pd,0,inf);
    r = random(t,length(v),1);
    tGsel.mN(ind)=mean(r);
    tGsel.sN(ind)=std(r);    
    [~,pN,ksstatN,~]= kstest(v,'CDF',t);
    tGsel.pN(ind)=pN;
    tGsel.ksstatN(ind)=ksstatN;    
    catch
        disp(ind)
        disp("Truncated Gaussian")
        disp("No convergence")
        disp("Fit Half Gaussian")
        disp(tGsel.np(ind))
        
        pd = fitdist(v','HalfNormal');
        tGsel.muN(ind)=pd.mu;
        tGsel.sigmaN(ind)=pd.sigma;
        r = random(pd,length(v),1);
        tGsel.mN(ind)=mean(r);
        tGsel.sN(ind)=std(r);    
        [~,pN,ksstatN,~]= kstest(v,'CDF',pd);
        tGsel.pN(ind)=pN;
        tGsel.ksstatN(ind)=ksstatN; 
        
    end    
    
    
     % exponential
    [muhat,muci] = expfit(v);
    tGsel.muE(ind)=muhat(1);
    pd = makedist('Exponential', muhat(1));
    r = random(pd,length(v),1);    
    tGsel.mE(ind)=mean(r);
    tGsel.sE(ind)=std(r);    
    [~,pE,ksstatE,~]= kstest(v,'CDF',pd);
    tGsel.pE(ind)=pE;
    tGsel.ksstatE(ind)=ksstatE;    
    
    
    
    
    % uniform
    [ahat,bhat] = unifit(v);
    tGsel.aU(ind)=ahat;
    tGsel.bU(ind)=bhat;
    pd=makedist('Uniform', ahat, bhat);
    r = random(pd,length(v),1);
    tGsel.mU(ind)=mean(r);
    tGsel.sU(ind)=std(r);    
    [~,pU,ksstatU,~]= kstest(v,'CDF',pd);
    tGsel.pU(ind)=pU;
    tGsel.ksstatU(ind)=ksstatU;    
    
end
