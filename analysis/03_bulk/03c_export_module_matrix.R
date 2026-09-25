
load("<ANALYSIS_ROOT>/fuxian_test/day2_1_bulk/rds/bulk_GSE132903_day2_1_objects.Rdata")
rownames(expr_gene) <- toupper(rownames(expr_gene))
genes <- strsplit("RHBDF2,SLC38A11,NOTCH3,PDGFRB,PTH1R,DGKB,GPC5,RNF152,GRM8,SLC6A12,INPP4B,CPM,CNTN4,GJC1,PDE7B,PLCE1,GRM3,CLMN,COL4A3,CACNA1C", ",")[[1]]
present <- intersect(genes, rownames(expr_gene))
mat <- t(expr_gene[present,,drop=FALSE])
df <- data.frame(sample_id=rownames(mat), group=as.character(group), mat, check.names=FALSE)
write.csv(df, "<ANALYSIS_ROOT>/fuxian_test/day6_GEM_method_enhancement_screen/tmp/bulk_GEM_gene_matrix.csv", row.names=FALSE, fileEncoding="UTF-8")
