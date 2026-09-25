.libPaths(c('<RLIB>','<R_SITE_LIBRARY>',.libPaths()))
suppressPackageStartupMessages({library(qs);library(miloR);library(Seurat);library(Matrix)})
root<-'<ANALYSIS_ROOT>/fuxian_test'
out<-'<WORKDIR>/results'
for(n in c('Cerebrovascular','Astrocytes','Microglia')) {
 p<-file.path(root,'day2_4_5_6_NVU_milor/rds',paste0('ckpt_',n,'_ADCN_Nature_v2_milor_final.qs'))
 x<-qread(p)
 write.csv(as.matrix(nhoodCounts(x$milo)),file.path(out,paste0(n,'_nhood_counts.csv')))
 write.csv(x$da,file.path(out,paste0(n,'_original_nhood_results.csv')),row.names=FALSE)
 write.csv(x$nhood_meta,file.path(out,paste0(n,'_nhood_metadata.csv')),row.names=FALSE)
 write.csv(x$sample_design,file.path(out,paste0(n,'_original_design.csv')),row.names=FALSE)
 print(x$parameters); print(head(x$nhood_meta,2)); print(table(x$da$Subtype))
 rm(x);gc()
}
x<-qread(file.path(root,'day2_2_subtype/rds/scRNA_cereb_annotated.qs'))
print(Assays(x));print(colnames(x@meta.data));print(dim(x));print(Layers(x[['RNA']]))
a<-if('decontX'%in%Assays(x)) 'decontX' else 'RNA'
m<-LayerData(x,assay=a,layer='counts');md<-x@meta.data
id<-factor(md$SampleID);if(length(id)==0)id<-factor(md$orig.ident)
mm<-sparse.model.matrix(~0+id);colnames(mm)<-levels(id)
pb<-m%*%mm
write.csv(as.matrix(pb),file.path(out,'vascular_pseudobulk_counts.csv'))
write.csv(data.frame(sample_id=levels(id),n_cells=as.integer(table(id))),file.path(out,'vascular_pseudobulk_ncells.csv'),row.names=FALSE)
write.csv(md,file.path(out,'vascular_cell_metadata.csv'))
writeLines(c(paste('assay',a),paste('genes',nrow(m)),paste('nuclei',ncol(m))),file.path(out,'vascular_pseudobulk_assay.txt'))
