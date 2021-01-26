function [tGsel]=groupstatistics_inchikey(db, tGsel)
tGsel.meanlog=zeros(height(tGsel),1);
tGsel.stdlog=zeros(height(tGsel),1);
tGsel.mean=zeros(height(tGsel),1);
tGsel.std=zeros(height(tGsel),1);
tGsel.numorganisms=zeros(height(tGsel),1);
tGsel.organisms=cell(height(tGsel),1);

for ind=1:height(tGsel)
    eclabel=tGsel.EC(ind);
    InChIKeylabel=tGsel.InChIKey(ind);
    
    filter=and(db.EC==eclabel,db.InChIKey==InChIKeylabel);
    tsel=db(filter,:);
    tGsel.meanlog(ind)=mean(log(tsel.KM));
    tGsel.stdlog(ind)=std(log(tsel.KM));
    tGsel.mean(ind)=mean(tsel.KM);
    tGsel.std(ind)=std(tsel.KM);
    tGsel.numorganisms(ind)=length(unique(tsel.organism));
    tGsel.organisms{ind}=strjoin(string(cellstr(unique(tsel.organism))'), '/');
end