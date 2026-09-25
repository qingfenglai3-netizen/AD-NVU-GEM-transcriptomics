.libPaths(c('<R_SITE_LIBRARY>','<R_USER_LIBRARY>',.libPaths()))
suppressPackageStartupMessages({library(Biobase);library(limma)})
out<-'<WORKDIR>/results'
load('<ANALYSIS_ROOT>/bulk_reference/GSE132903_eSet.Rdata')
load('<ANALYSIS_ROOT>/bulk_reference/GPL10558_bioc.rda')
original<-new.env();load('<ANALYSIS_ROOT>/fuxian_test/day2_1_bulk/rds/bulk_GSE132903_day2_1_objects.Rdata',envir=original)
tr<-read.csv(file.path(out,'discovery_GSE132903_features.csv'),check.names=FALSE);genes<-colnames(tr)[-(1:2)]
x<-exprs(gset[[1]]);a<-GPL10558_bioc;colnames(a)<-c('probe_id','symbol');a$symbol<-trimws(as.character(a$symbol))
a<-a[a$probe_id %in% rownames(x)& !is.na(a$symbol)& a$symbol!='',];a<-a[a$symbol %in% genes,]
stopifnot(!anyDuplicated(a$probe_id),!anyDuplicated(tr$sample_id),setequal(colnames(x),tr$sample_id))
y<-factor(tr$group,levels=c('CN','AD'));x<-x[,tr$sample_id,drop=FALSE]
fit<-eBayes(lmFit(x[a$probe_id,,drop=FALSE],model.matrix(~y)))
effects<-topTable(fit,coef=2,number=Inf,sort.by='none');a$logFC_all195<-effects[a$probe_id,'logFC']
a$archived_selected<-a$probe_id %in% original$all_diff$probe_id
write.csv(a,file.path(out,'ML_probe_lineage.csv'),row.names=FALSE)
avg<-t(sapply(genes,function(g) colMeans(x[a$probe_id[a$symbol==g],,drop=FALSE])))
result<-data.frame(sample_id=tr$sample_id,group=tr$group,t(avg),check.names=FALSE)
write.csv(result,file.path(out,'discovery_GSE132903_labelblind_mean_features.csv'),row.names=FALSE)
summary<-do.call(rbind,lapply(genes,function(g){aa<-a[a$symbol==g,];v<-x[aa$probe_id[aa$archived_selected],];data.frame(gene=g,n_probes=nrow(aa),archived_probe=paste(aa$probe_id[aa$archived_selected],collapse=';'),max_abs_input_difference=max(abs(v-tr[[g]])),mean_aggregation_difference=max(abs(result[[g]]-tr[[g]])),selected_is_max_abs_effect=all(abs(aa$logFC_all195[aa$archived_selected])>=max(abs(aa$logFC_all195))-1e-10))}))
stopifnot(all(summary$max_abs_input_difference<1e-10),all(summary$selected_is_max_abs_effect))
write.csv(summary,file.path(out,'ML_probe_lineage_summary.csv'),row.names=FALSE)
print(summary);print(sessionInfo())
