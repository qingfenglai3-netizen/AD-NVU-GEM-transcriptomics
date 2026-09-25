suppressPackageStartupMessages({library(edgeR);library(limma)})
set.seed(20260910)
base<-normalizePath(Sys.getenv('GEM_REPO',getwd()),winslash='/');input<-file.path(base,'data');out<-Sys.getenv('GEM_OUTPUT',file.path(base,'results'));dir.create(out,recursive=TRUE,showWarnings=FALSE)
dm<-read.csv(file.path(input,'donor_APOE_manifest.csv'));nc<-read.csv(file.path(input,'vascular_pseudobulk_ncells.csv'))
dm<-merge(dm,nc,by.x='SampleID',by.y='sample_id');dm<-dm[dm$Dataset=='D1'&dm$n_cells>=20,];rownames(dm)<-dm$SampleID;dm$Group<-factor(dm$Group,levels=c('CN','AD'));dm$e4carrier<-factor(ifelse(dm$APOE4_dosage>0,'carrier','noncarrier'),levels=c('noncarrier','carrier'))
write.csv(as.data.frame(table(dm$e4carrier,dm$Group)),file.path(out,'APOE4_by_diagnosis_complete_case.csv'),row.names=FALSE)
pb<-as.matrix(read.csv(file.path(input,'vascular_pseudobulk_counts.csv'),row.names=1,check.names=FALSE));genes<-read.csv(file.path(input,'vascular_DE_cohort_adjusted.csv'))$gene
ys<-DGEList(pb[genes,rownames(dm)]);ys<-calcNormFactors(ys)
specs<-list(diagnosis_only=~Group,sex_adjusted=~sex+Group,sex_APOE_dosage=~sex+APOE2_dosage+APOE4_dosage+Group,sex_APOE4_carrier=~sex+e4carrier+Group,interaction=~sex+Group*e4carrier,noncarrier=~sex+Group,carrier=~sex+Group)
rows<-list();coefs<-list()
for(model in names(specs)){
 dd<-if(model%in%c('noncarrier','carrier'))dm[dm$e4carrier==model,] else dm
 des<-model.matrix(specs[[model]],dd);stopifnot(qr(des)$rank==ncol(des));y<-ys[,rownames(dd),keep.lib.sizes=TRUE]
 y<-estimateDisp(y,des,robust=TRUE);fit<-glmQLFit(y,des,robust=TRUE);term<-if(model=='interaction')'GroupAD:e4carriercarrier' else 'GroupAD';test<-glmQLFTest(fit,coef=term);tab<-topTags(test,n=Inf)$table;g<-match('GEM',rownames(fit));tt<-tab['GEM',]
 # Approximate quasi-likelihood Wald interval for the unshrunk coefficient.
 if(model=='diagnosis_only'){print(names(fit));print(str(fit$s2.post));print(str(fit$var.post))}
 mu<-fit$fitted.values[g,];phi<-if(length(fit$dispersion)==1)fit$dispersion else fit$dispersion[g];w<-mu/(1+phi*mu);post<-if(!is.null(fit$var.post))fit$var.post else fit$s2.post;pv<-if(length(post)==1)post else post[g];cov<-solve(crossprod(des,des*w))*pv;se<-sqrt(cov[term,term])/log(2)
 beta<-fit$unshrunk.coefficients[g,term]/log(2);df<-test$df.total[g];if(length(df)==0 || !is.finite(df))df<-max(1,nrow(des)-ncol(des))
 ci<-beta+c(-1,1)*qt(.975,df)*se
 rows[[length(rows)+1]]<-data.frame(model=model,term=term,n_AD=sum(dd$Group=='AD'),n_CN=sum(dd$Group=='CN'),edgeR_logFC=tt$logFC,unshrunk_logFC=beta,approx_QL_SE=se,approx_QL_CI_low=ci[1],approx_QL_CI_high=ci[2],PValue=tt$PValue,FDR=tt$FDR,design_rank=qr(des)$rank,design_columns=ncol(des),condition_number=kappa(des),df_approx=df)
 write.csv(cbind(sample_id=rownames(dd),des),file.path(out,paste0('APOE_nested_design_',model,'.csv')),row.names=FALSE)
}
ans<-do.call(rbind,rows);write.csv(ans,file.path(out,'APOE_nested_stratified_GEM.csv'),row.names=FALSE);print(ans)
# Donor influence (Cook's distance) on the GEM response.
z<-dm;z$GEM_logCPM<-as.numeric(cpm(ys,log=TRUE,prior.count=2)['GEM',]);m<-lm(GEM_logCPM~sex+APOE2_dosage+APOE4_dosage+Group,z)
inf<-data.frame(sample_id=rownames(z),GEM_logCPM=z$GEM_logCPM,Group=z$Group,n_cells=z$n_cells,Cooks_distance=cooks.distance(m),hat=hatvalues(m))
write.csv(inf,file.path(out,'APOE_GEM_donor_influence.csv'),row.names=FALSE)
write.csv(z,file.path(out,'APOE_GEM_complete_case_values.csv'),row.names=FALSE)
