function tlog=script_SK_nutrients(tlog, fooddatabase0910selclean, muf,sf)
vTrunc=100;
tlog.muf=muf;
tlog.sf=sf;
tlog.Tm=zeros(height(tlog),1);
tlog.sk=zeros(height(tlog),1);

tlog.TmL=zeros(height(tlog),1);
tlog.skL=zeros(height(tlog),1);

tlog.TmG=zeros(height(tlog),1);
tlog.skG=zeros(height(tlog),1);

tlog.TmW=zeros(height(tlog),1);
tlog.skW=zeros(height(tlog),1);

tlog.TmN=zeros(height(tlog),1);
tlog.skN=zeros(height(tlog),1);



for ind=1:size(fooddatabase0910selclean,1)
    v=fooddatabase0910selclean(ind,:);
    v=v(v>0);
    if or(or(or(tlog.muf(ind)==0, tlog.sf(ind)==0), isnan(tlog.muf(ind))), isnan(tlog.sf(ind)))
        disp(ind)
        disp("No content")
        tlog.Tm(ind)=nan;
        tlog.sk(ind)=nan;
        tlog.TmL(ind)=nan;
        tlog.skL(ind)=nan;
        tlog.TmG(ind)=nan;
        tlog.skG(ind)=nan; 
        tlog.TmW(ind)=nan;
        tlog.skW(ind)=nan; 
        tlog.TmN(ind)=nan;
        tlog.skN(ind)=nan;  
        
    else   
  % real data   
    tlog.Tm(ind)=moment(log(v),3);
    tlog.sk(ind)=skewness(log(v)); 
  % lognormal
    pdf_truncln = @(v,mu,sigma) lognpdf(v,mu,sigma) ./ logncdf(vTrunc,mu,sigma);
    [lognp,lognci] = lognfit(v);
    start = [lognp(1),lognp(2)];
    try
    [paramEstsL,paramCIsL] = mle(v, 'pdf',pdf_truncln, 'start',start, 'lower',[-Inf 0]);
    pd = makedist('Lognormal', paramEstsL(1), paramEstsL(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    tlog.TmL(ind)=moment(log(r),3);
    tlog.skL(ind)=skewness(log(r));    
    catch
        disp(ind)
        disp("Truncated Lognormal")
        disp("No convergence")
        tlog.TmL(ind)=nan;
        tlog.skL(ind)=nan;     
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
    tlog.TmG(ind)=moment(log(r),3);
    tlog.skG(ind)=skewness(log(r));    
    catch
        disp(ind)
        disp("Truncated Gamma")
        disp("No convergence")
        tlog.TmG(ind)=nan;
        tlog.skG(ind)=nan;    
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
    tlog.TmW(ind)=moment(log(r),3);
    tlog.skW(ind)=skewness(log(r));  

    catch
        disp(ind)
        disp("Truncated Weibull")
        disp("No convergence")
        tlog.TmW(ind)=nan;
        tlog.skW(ind)=nan;   
    end     
    
    % gaussian
    pdf_truncnorm = @(v,mu,sigma) normpdf(v,mu,sigma) ./ (normcdf(vTrunc,mu,sigma)-normcdf(0,mu,sigma));
    start = [mean(v),std(v)];
    try
    [paramEstsN,paramCIsN] = mle(v, 'pdf',pdf_truncnorm, 'start',start, 'lower',[0 0]);  
    pd = makedist('Normal', paramEstsN(1), paramEstsN(2));
    t = truncate(pd,0,vTrunc);
    r = random(t,length(v),1);
    tlog.TmN(ind)=moment(log(r),3);
    tlog.skN(ind)=skewness(log(r));  

    catch
        disp(ind)
        disp("Truncated Gaussian")
        disp("No convergence")
        disp("Fit Half Gaussian")

        pd = fitdist(v','HalfNormal');
        r = random(pd,length(v),1);
        tlog.TmN(ind)=moment(log(r),3);
        tlog.skN(ind)=skewness(log(r));  
        
    end
    end
end