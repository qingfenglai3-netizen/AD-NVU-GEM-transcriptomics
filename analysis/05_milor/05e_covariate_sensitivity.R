.libPaths(c('<RLIB>','<R_SITE_LIBRARY>',.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma);library(qs);library(miloR)})
set.seed(20260910)
root<-'<ANALYSIS_ROOT>/fuxian_test'
out<-'<WORKDIR>/results'
dm<-read.csv(file.path(out,'donor_APOE_manifest.csv'))
summary<-list()
for(n in c('Cerebrovascular','Astrocytes','Microglia')) {
 x<-qread(file.path(root,'day2_4_5_6_NVU_milor/rds',paste0('ckpt_',n,'_ADCN_Nature_v2_milor_final.qs')))
 sd<-dm[match(colnames(nhoodCounts(x$milo)),dm$SampleID),];rownames(sd)<-sd$SampleID
 sd$Group<-factor(sd$Group,levels=c('CN','AD'));sd$Dataset<-factor(sd$Dataset)
 for(model in c('diagnosis_only','cohort_adjusted','D1_APOE_sex_adjusted')){
  if(model=='D1_APOE_sex_adjusted') {
   # Same neighbourhood graph; D1 donor counts tested in edgeR with TMM normalization.
   dd<-sd[sd$Dataset=='D1',];f<-~APOE2_dosage+APOE4_dosage+sex+Group
  } else {dd<-sd;f<-if(model=='diagnosis_only') ~Group else ~Dataset+Group}
  des<-model.matrix(f,dd);stopifnot(qr(des)$rank==ncol(des))
  y<-DGEList(as.matrix(nhoodCounts(x$milo))[,rownames(dd),drop=FALSE]);y<-calcNormFactors(y,method='TMM')
  y<-estimateDisp(y,des);fit<-glmQLFit(y,des,robust=TRUE,legacy=TRUE);test<-glmQLFTest(fit,coef='GroupAD')
  tab<-topTags(test,n=Inf,sort.by='none')$table;tab$Nhood<-x$da$Nhood;tab$Subtype<-x$da$NhoodSubtype
  tab$FDR<-p.adjust(tab$PValue,'BH');tab$model<-model;tab$compartment<-n
  write.csv(tab,file.path(out,paste0(n,'_MiloR_',model,'.csv')),row.names=FALSE)
  write.csv(cbind(sample_id=rownames(dd),des),file.path(out,paste0(n,'_design_',model,'.csv')),row.names=FALSE)
  summary[[length(summary)+1]]<-data.frame(compartment=n,model=model,n_donors=nrow(dd),n_AD=sum(dd$Group=='AD'),n_CN=sum(dd$Group=='CN'),n_tests=nrow(tab),n_FDR10=sum(tab$FDR<.1),min_FDR=min(tab$FDR),max_abs_diff_original_logFC=if(model=='diagnosis_only')max(abs(tab$logFC-x$da$logFC)) else NA)
 }
 rm(x);gc()
}
write.csv(do.call(rbind,summary),file.path(out,'MiloR_covariate_summary.csv'),row.names=FALSE)
pb<-as.matrix(read.csv(file.path(out,'vascular_pseudobulk_counts.csv'),row.names=1,check.names=FALSE))
nc<-read.csv(file.path(out,'vascular_pseudobulk_ncells.csv'));keep<-nc$sample_id[nc$n_cells>=20]
pb<-pb[,keep,drop=FALSE];sd<-dm[match(keep,dm$SampleID),];rownames(sd)<-sd$SampleID
sd$Group<-factor(sd$Group,levels=c('CN','AD'));sd$Dataset<-factor(sd$Dataset)
res<-list()
frozen_genes<-NULL
for(model in c('cohort_adjusted','D1_diagnosis_only','D1_APOE_sex_adjusted')){
 if(grepl('^D1',model)){dd<-sd[sd$Dataset=='D1',];f<-if(model=='D1_diagnosis_only')~Group else ~APOE2_dosage+APOE4_dosage+sex+Group}else{dd<-sd;f<-~Dataset+Group}
 des<-model.matrix(f,dd);stopifnot(qr(des)$rank==ncol(des))
 y<-DGEList(pb[,rownames(dd),drop=FALSE]);kg<-if(model=='cohort_adjusted')filterByExpr(y,des,min.count=10) else rownames(y)%in%frozen_genes;y<-y[kg,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,des,robust=TRUE);fit<-glmQLFit(y,des,robust=TRUE);test<-glmQLFTest(fit,coef='GroupAD');tab<-topTags(test,n=Inf)$table;tab$gene<-rownames(tab)
 if(model=='cohort_adjusted')frozen_genes<-tab$gene
 write.csv(tab,file.path(out,paste0('vascular_DE_',model,'.csv')),row.names=FALSE)
 gem<-tab[tab$gene=='GEM',];gem$model<-model;gem$n_AD<-sum(dd$Group=='AD');gem$n_CN<-sum(dd$Group=='CN');gem$n_genes<-nrow(tab);gem$n_FDR05<-sum(tab$FDR<.05);res[[length(res)+1]]<-gem
 write.csv(cbind(sample_id=rownames(dd),des),file.path(out,paste0('vascular_DE_design_',model,'.csv')),row.names=FALSE)
 write.csv(cpm(y,log=TRUE,prior.count=2)['GEM',,drop=FALSE],file.path(out,paste0('vascular_GEM_logCPM_',model,'.csv')))
}
write.csv(do.call(rbind,res),file.path(out,'GEM_APOE_sensitivity.csv'),row.names=FALSE)
print(do.call(rbind,summary));print(do.call(rbind,res));capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
