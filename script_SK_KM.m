function tGsel=script_SK_KM(tGsel, db)


tGsel.sk=zeros(height(tGsel),1);
tGsel.skL=zeros(height(tGsel),1);
tGsel.skG=zeros(height(tGsel),1);
tGsel.skW=zeros(height(tGsel),1);
tGsel.skN=zeros(height(tGsel),1);
tGsel.skE=zeros(height(tGsel),1);
tGsel.skU=zeros(height(tGsel),1);


for ind=1:height(tGsel)
    eclabel=tGsel.EC(ind);
    sublabel=tGsel.substrate(ind);
    
    filter=and(db.EC==eclabel,db.substrate==sublabel);
    tsel=db(filter,:);
    v=tsel.KM;
    tGsel.sk(ind)=skewness(log(v)); 
    
     % lognormal
    [lognp,lognci] = lognfit(v);
    pd = makedist('Lognormal', lognp(1), lognp(2));
    r = random(pd,length(v),1);
    tGsel.skL(ind)=skewness(log(r));  

        
    
    % gamma
    [gammap,gammaci] = gamfit(v);  
    pd = makedist('Gamma', gammap(1), gammap(2));
    r = random(pd,length(v),1);
    tGsel.skG(ind)=skewness(log(r)); 

    
 
    % weibull
    [wblp,wblci] = wblfit(v);
    pd = makedist('Weibull', wblp(1), wblp(2));
    r = random(pd,length(v),1);
    tGsel.skW(ind)=skewness(log(r));  
    
    
    %gaussian   
    pdf_truncnorm = @(v,mu,sigma) normpdf(v,mu,sigma) ./ (normcdf(inf,mu,sigma)-normcdf(0,mu,sigma));
    start = [mean(v),std(v)];
    try
    [paramEstsN,paramCIsN] = mle(v, 'pdf',pdf_truncnorm, 'start',start, 'lower',[0 0]);   
    pd = makedist('Normal', paramEstsN(1), paramEstsN(2));
    t = truncate(pd,0,inf);
    r = random(t,length(v),1);    
    catch
        disp(ind)
        disp("Truncated Gaussian")
        disp("No convergence")
        disp("Fit Half Gaussian")
        disp(tGsel.np(ind))
        
        pd = fitdist(v','HalfNormal');
        r = random(pd,length(v),1);
        
    end    
    tGsel.skN(ind)=skewness(log(r));
    
     % exponential
    [muhat,muci] = expfit(v);
    pd = makedist('Exponential', muhat(1));
    r = random(pd,length(v),1);
    tGsel.skE(ind)=skewness(log(r));
 
    
    % uniform
    [ahat,bhat] = unifit(v);
    pd=makedist('Uniform', ahat, bhat);
    r = random(pd,length(v),1);
    tGsel.skU(ind)=skewness(log(r));
 
end