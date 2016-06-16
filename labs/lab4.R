comp.fun<-function(grm) {#grm needs to be a data.frame with first column being the first id, secnod column being the second id (with names id1 and id2) and the third column the genetic relatedness value (name=grm)
lf<-c(ceu="CEU.Illu550K.list",mes="MESTIZOS.Illu550K.list",nat="NATIVES.Illu550K.list",yri="YRI.Illu550K.list")
L<-list()
for (i in 1:length(lf)) {
    read.table(paste("~/mexdat/",lf[i],sep=""),stringsAsFactors=FALSE,colClasses="character")->id
    id[,1:2]->id
    names(id)<-c("fam","id")
    names(lf)[i]->id$grp
    gsub("_sex","",id$id)->id$id
    id->L[[i]]
}
data.frame(do.call("rbind",L))->grp
NULL->grp$fam #don't seem to need this
#
names(grp)->nms
names(grp)<-paste(nms,"1",sep="")
merge(grp,grm,all.y=TRUE)->x
names(grp)<-paste(nms,"2",sep="")
merge(grp,x,all.y=TRUE)->x
nms->names(grp)
#
x[x$id1!=x$id2,]->x2
unique(grp$grp)->grps
mat1<-mat2<-matrix(NA,length(grps),length(grps))
for (i in 1:length(grps)) for (j in 1:length(grps)) {
    x2[x2$grp1==grps[i] & x2$grp2==grps[j],]->tmp
    mean(tmp$grm,na.rm=TRUE)->mat1[i,j]
    sd(tmp$grm,na.rm=TRUE)->mat2[i,j]
}
grps->rownames(mat1)->colnames(mat1)->rownames(mat2)->colnames(mat2)
list(x2,mat1,mat2)
}


Lab #4
For this labe, we are going to focus on the Mexico data (doi: 10.1073/pnas.0903045106). Recall that this data has both Mestizos and Amerindians, so we'll want to track that. 
Let's first prune for LD so that things are a little simpler.
cd ~/wd/
plink --bfile ~/mexdat/CEU.YRI.MESTIZOS.NATIVO.550K.Common.1-22.v2.flipped --indep-pairwise 50 5 0.2 --out mx_le #why are correlations being computed within chrosomome?
plink --bfile ~/mexdat/CEU.YRI.MESTIZOS.NATIVO.550K.Common.1-22.v2.flipped --extract mx_le.prune.in --make-bed --out mx_le


A. GCTA/plink GRM
1. Use the "--make-grm-gz" option in plink to make a genetic relatedness matrix (GRM). 
plink --bfile mx_le --make-grm-gz  --out mx_le

2. Take a look at how the resulting values are distributed. How many are there? Can you make sense of that number?
gunzip -c mx_le.grm.gz > mx_le.grm

#R
read.table("mx_le.grm",stringsAsFactors=FALSE)->grm
read.table("mx_le.grm.id",stringsAsFactors=FALSE)->id
NULL->grm[,3]
names(grm)<-c("id1","id2","grm")
id[,2][grm[,1] ] -> grm[,1]
id[,2][grm[,2] ] -> grm[,2]
#id[,1][grm[,1] ] -> grm$fam1
#id[,1][grm[,2] ] -> grm$fam2
comp.fun(grm)->gcta.comp


#hist(grm$grm,breaks=100) #this is a problem!


B. plink pi-hat
plink --bfile mx_le --genome --out mx_le

#R
read.table("mx_le.genome",header=TRUE)->gen
summary(gen$PI_HAT)
#hist(gen$PI_HAT,breaks=100)
gen[,c("IID1","IID2","DST")]->tmp
names(tmp)<-c("id1","id2","grm")
comp.fun(tmp)->plink.comp

C. king
king -b mx_le.bed --kinship    
awk '{print $2" "$4" "$8}' king.kin0 > king_reduced.kin0 #what is in this file?

read.table("king_reduced.kin0",header=TRUE)->k
names(k)<-c("id1","id2","grm")
comp.fun(k)->kin.comp

D. reap
#admixture mx_le.bed  --cv 4 &
admixture mx_le.bed   3 

#need to combine IDs and admixture proportions
read.table("mx_le.3.Q")->ind_pr
read.table("mx_le.fam")->fam
cbind(fam[,1:2],ind_pr)->tmp
write.table(tmp,file="admix.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

plink --bfile mx_le --recode transpose 12 --out mx_leT
REAP -g mx_leT.tped -p mx_leT.tfam -a admix.txt -f mx_le.3.P -k 3 -t -1 -m

awk '{print $2, $4, $9}' REAP_pairs_relatedness.txt > reap_pairs.txt
read.table("reap_pairs.txt",header=TRUE)->k
names(k)<-c("id1","id2","grm")
comp.fun(k)->reap.comp

#write some code to compare all
L<-list(gcta=gcta.comp[[1]],plink=plink.comp[[1]],king=kin.comp[[1]],reap=reap.comp[[1]])
for (i in 1:length(L)) {
    L[[i]]->zz
    nm<-names(L)[i]
    nm->names(zz)[5]
    fun<-function(x) paste(sort(x),collapse=".")
    apply(zz[,c("id1","id2")],1,fun)->zz$id
    apply(zz[,c("grp1","grp2")],1,fun)->zz$grp
    zz[,c("id","grp",nm)]->zz
    zz->L[[i]]
}
L[[1]]->df
for (i in 2:length(L)) merge(df,L[[i]],all=TRUE)->df

par(mfrow=c(2,2))
for (nm in c("gcta","plink","king","reap")) {
    df[[nm]]->z
    plot(density(z,na.rm=TRUE))
}
plot(df[,c("gcta","plink","king","reap")])

split(df,df$grp)->tmp
fun<-function(x) cor(x[,c("gcta","plink","king","reap")],use='p')
fun(df)
lapply(tmp,fun)

