function tlog=script_SD_nutrients(tlog, fooddatabase0910selclean, muf,sf)
vTrunc=100;
tlog.muf=muf;
tlog.sf=sf;
tlog.sd=zeros(height(tlog),1);
tlog.sdL=zeros(height(tlog),1);
tlog.sdG=zeros(height(tlog),1);
tlog.sdW=zeros(height(tlog),1);
tlog.sdN=zeros(height(tlog),1);

for ind=1:size(fooddatabase0910selclean,1)
    v=fooddatabase0910selclean(ind,:);
    v=v(v>0);
    if or(or(or(tlog.muf(ind)==0, tlog.sf(ind)==0), isnan(tlog.muf(ind))), isnan(tlog.sf(ind)))
    
        disp(ind)
        disp("No content")
        tlog.sd(ind)=nan;
        tlog.sdL(ind)=nan;
        tlog.sdG(ind)=nan; 
        tlog.sdW(ind)=nan;  
        tlog.sdN(ind)=nan;  
        
    else   
  % real data   
    tlog.sd(ind)=std(log(v)); 
  % lognormal
    pdf_truncln = @(v,mu,sigma) lognpdf(v,mu,sigma) ./ logncdf(vTrunc,mu,sigma);
    [lognp,lognci] = lognfit(v);
    start = [lognp(1),lognp(2)];
    try
    [paramEstsL,paramCIsL] = mle(v, 'pdf',pdf_truncln, 'start',start, 'lower',[-Inf 0]);
    pd = makedist('Lognormal', paramEstsL(1), paramEstsL(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    tlog.sdL(ind)=std(log(r));    
    catch
        disp(ind)
        disp("Truncated Lognormal")
        disp("No convergence")
        tlog.sdL(ind)=nan;     
    end     
    
    
 % gamma   
     pdf_truncgamma = @(v,a,b) gampdf(v,a,b) ./ gamcdf(vTrunc,a,b);
    [gammap,gammaci] = gamfit(v);
    start = [gammap(1),gammap(2)];   
    try
    [paramEstsG,paramCIsG] = mle(v, 'pdf',pdf_truncgamma, 'start',start, 'lower',[0 0]);
     pd = makedist('Gamma', paramEstsG(1), paramEstsG(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    tlog.sdG(ind)=std(log(r));    
    catch
        disp(ind)
        disp("Truncated Gamma")
        disp("No convergence")
        tlog.sdG(ind)=nan;    
    end   
    
    % weibull
    pdf_truncwbl = @(v,a,b) wblpdf(v,a,b) ./ wblcdf(vTrunc,a,b);
    [wblp,wblci] = wblfit(v);
    start = [wblp(1),wblp(2)];
    try
    [paramEstsW,paramCIsW] = mle(v, 'pdf',pdf_truncwbl, 'start',start, 'lower',[0 0]);
    pd = makedist('Weibull', paramEstsW(1), paramEstsW(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    tlog.sdW(ind)=std(log(r));  

    catch
        disp(ind)
        disp("Truncated Weibull")
        disp("No convergence")
        tlog.sdW(ind)=nan;   
    end     
    
    % gaussian
    pdf_truncnorm = @(v,mu,sigma) normpdf(v,mu,sigma) ./ (normcdf(vTrunc,mu,sigma)-normcdf(0,mu,sigma));
    start = [mean(v),std(v)];
    try
    [paramEstsN,paramCIsN] = mle(v, 'pdf',pdf_truncnorm, 'start',start, 'lower',[0 0]);  
    pd = makedist('Normal', paramEstsN(1), paramEstsN(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    tlog.sdN(ind)=std(log(r));  

    catch
        disp(ind)
        disp("Truncated Gaussian")
        disp("No convergence")
        disp("Fit Half Gaussian")

        pd = fitdist(v','HalfNormal');
        r = random(pd,length(v),1);
        tlog.sdN(ind)=std(log(r));  
        
    end
    end
end