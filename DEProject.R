#microarray data: GSE110978 GSE117763 GSE137176
#from GSE110978: 8 mice , 4 young, 4 old -> fibroblasts
#from GSE117763: 8 mice, 4 young, 4 old -> fibroblasts
#from GSE137176: 6 mice, 3 young, 3 old -> epidermal 
#goal is to find DE genes in epidermal or fibroblasts between two conditions: young and old 

library(tidyverse)
library(affy)
library(limma)
library(Glimma)

#Step 1: Read in the three folders of cel files and convert to ExpressionSetObjects. 
setwd('/Users/Emix/Desktop/DESeq/GSE110978_RAW')
study1 <- ReadAffy()
fibroblast1 <- rma(study1)

setwd('/Users/Emix/Desktop/DESeq/GSE117763_RAW')
study2 <- ReadAffy()
fibroblast2 <- rma(study2)

setwd('/Users/Emix/Desktop/DESeq/GSE137176_RAW')
study3 <- ReadAffy()
epidermal <- rma(study3)


#Step 2: Add the condition column to each pData and explore what data looks like. 
pData(fibroblast1)$condition <- c(rep('young', 4), rep('old',4))
pdf('density_fb1.pdf')
plotDensities(fibroblast1, group = pData(fibroblast1)[, "condition"], legend = FALSE)
dev.off()

pData(fibroblast2)$condition <- c(rep(c('young', 'old'), 4))
pdf('density_fb2.pdf')
plotDensities(fibroblast2, group = pData(fibroblast2)[, "condition"], legend = FALSE)
dev.off()

pData(epidermal)$condition <- c(rep('young', 3), rep('old',3))
pdf('density_epi.pdf')
plotDensities(epidermal, group = pData(epidermal)[, "condition"], legend = FALSE)
dev.off()

#density plots are not very right-skewed therefore many genes have expression!

#Step 3: Create design matrix for each. 
fb1_design <- model.matrix(~0 + condition, data = pData(fibroblast1))
fb2_design <- model.matrix(~0 + condition, data = pData(fibroblast2))
epi_design <- model.matrix(~0 + condition, data = pData(epidermal))
#check design matrices 
colSums(fb1_design) #4 young, 4 old
colSums(fb2_design) #4 young, 4 old
colSums(epi_design)#3 young, 3 old

#Step 4: Create contrast matrix for each. 
fb1_cm <- makeContrasts(status = conditionold - conditionyoung, 
                    levels = fb1_design)

fb2_cm <- makeContrasts(status = conditionold - conditionyoung, 
                        levels = fb2_design)

epi_cm <- makeContrasts(status = conditionold - conditionyoung, 
                       levels = epi_design)

#Step 5: Perform the rest of limma pipeline and view results. 
fb_fit1 <-lmFit(fibroblast1, fb1_design)
fb_fit1 <- contrasts.fit(fb_fit1, fb1_cm)
fb_fit1 <- eBayes(fb_fit1) 
fb1_results <- decideTests(fb_fit1) #summarize results 
summary(fb1_results) #2481 down, 3435 up regulated genes 

fb_fit2 <-lmFit(fibroblast2, fb2_design)
fb_fit2 <- contrasts.fit(fb_fit2, fb2_cm)
fb_fit2 <- eBayes(fb_fit2) 
fb2_results <- decideTests(fb_fit2) #summarize results 
summary(fb2_results) #108 down, 120 up regulated genes 


epi_fit <-lmFit(epidermal, epi_design)
epi_fit <- contrasts.fit(epi_fit, epi_cm)
epi_fit <- eBayes(epi_fit) 
epi_results <- decideTests(epi_fit) #summarize results 
summary(epi_results) #4001 down, 2175 up regulated genes 


#Step 6: Pull gene IDs and gene symbols for results and visualize them. 
annotation(fibroblast1)
annotation(fibroblast2)
annotation(epidermal)

library(htmg430pm.db)
library(mogene10sttranscriptcluster.db)
library(annotate)
library(EnhancedVolcano)

fb1_annotations<- select(htmg430pm.db, keys = rownames(fb_fit1), columns = c("ENTREZID", "SYMBOL"))
fb1_annotations <- fb1_annotations[!duplicated(fb1_annotations$PROBEID), ]
#one probe ID mapped to multiple entrezIDs // removing entries after first occurrence to simplify analysis 
stats_fb1 <- topTable(fb_fit1, coef= 'status', number = nrow(fb_fit1), sort.by = 'none')
stats_fb1$SYMBOL <- fb1_annotations$SYMBOL
  

#set parameters so can see the v significant genes 
pdf('fb1volcano.pdf')
EnhancedVolcano(stats_fb1, lab = stats_fb1$SYMBOL, x = 'logFC', y='P.Value', 
                pCutoff = 0.05, FCcutoff = 1)

EnhancedVolcano(test, lab = test$SYMBOL, x = 'logFC', y='P.Value', 
                pCutoff = 0.05, FCcutoff = 1)
dev.off()

#histogram to check distribution 
pdf('hist_fb1.pdf')
hist(stats_fb1[, "P.Value"])
dev.off()


fb2_annotations<- select(htmg430pm.db, keys=rownames(fb_fit2), columns=c("ENTREZID", "SYMBOL")) 
fb2_annotations <- fb2_annotations[!duplicated(fb2_annotations$PROBEID), ]
stats_fb2 <- topTable(fb_fit2, coef= 'status', number = nrow(fb_fit2), sort.by = 'none')
stats_fb2$SYMBOL <- fb2_annotations$SYMBOL

pdf('fb2volcano.pdf')
EnhancedVolcano(stats_fb2, lab = stats_fb2$SYMBOL, x = 'logFC', y='P.Value', 
                pCutoff = 0.05, FCcutoff = 1)
dev.off()

pdf('hist_fb2.pdf')
hist(stats_fb2[, "P.Value"])
dev.off()

epi_annotations <- select(mogene10sttranscriptcluster.db, keys=rownames(epi_fit), columns=c("ENTREZID", "SYMBOL")) 
epi_annotations <- epi_annotations[!duplicated(epi_annotations$PROBEID), ]
stats_epi <- topTable(epi_fit, coef= 'status', number = nrow(epi_fit), sort.by = 'none')
stats_epi$SYMBOL <- epi_annotations$SYMBOL

pdf('epivolcano.pdf')
EnhancedVolcano(stats_epi, lab = stats_epi$SYMBOL, x = 'logFC', y='P.Value', 
                pCutoff = 0.05, FCcutoff = 1)
dev.off()

pdf('hist_epi.pdf')
hist(stats_epi[, "P.Value"])
dev.off()

#Step 7: Visualize differential gene expression results further and generate 
#HTML reports of the mean-difference plots for each dataset. 

glMDPlot(fb_fit1, counts = exprs(fibroblast1), anno = fb1_annotations, 
         status = fb1_results, groups = pData(fibroblast1)$condition, 
         main = 'Fibroblast DE Genes 1', 
         html = 'MD-Plot-Fib1')

glMDPlot(fb_fit2, counts = exprs(fibroblast2), anno = fb2_annotations, 
         status = fb2_results, groups = pData(fibroblast2)$condition, 
         main = 'Fibroblast DE Genes 2', 
         html = 'MD-Plot-Fib2')

glMDPlot(epi_fit, counts = exprs(epidermal), anno = epi_annotations, 
         status = epi_results, groups = pData(epidermal)$condition, 
         main = 'Epidermal DE Genes', 
         html = 'MD-Plot-Epi')


#Step 8: Link to KEGG and GO to obtain pathway information. 

enrich_fib1 <- kegga(fb_fit1, geneid = fb1_annotations$ENTREZID, species = 'Mm')
top_keggfib1 <- topKEGG(enrich_fib1, number = 5)
enrich_go_f1 <- goana(fb_fit1, geneid = fb1_annotations$ENTREZID, species = "Mm")
top_gofib1 <- topGO(enrich_go_f1, ontology = 'BP', number = 5)

write.csv(top_keggfib1, 'top_keggfib1.csv')
write.csv(top_gofib1, 'top_gofib1.csv')

enrich_fib2 <- kegga(fb_fit2, geneid = fb2_annotations$ENTREZID, species = 'Mm')
top_keggfib2<- topKEGG(enrich_fib2, number = 5)
enrich_go_f2 <- goana(fb_fit2, geneid = fb2_annotations$ENTREZID, species = "Mm")
top_gofib2 <- topGO(enrich_go_f2, ontology = 'BP', number = 5)

write.csv(top_keggfib2, 'top_keggfib2.csv')
write.csv(top_gofib2, 'top_gofib2.csv')

enrich_epi <- kegga(epi_fit, geneid = epi_annotations$ENTREZID, species = 'Mm')
top_keggepi<- topKEGG(enrich_epi, number = 5)
enrich_go_epi <- goana(epi_fit, geneid = epi_annotations$ENTREZID, species = "Mm")
top_goepi <- topGO(enrich_go_epi, ontology = 'BP', number = 5)

write.csv(top_keggepi, 'top_kegg_epi.csv')
write.csv(top_goepi, 'top_go_epi.csv')


#Step 9: Obtain list of top 10 distinct DE genes for each dataset 
#and visualize these results. 

top10_fib1 <- stats_fb1 %>% 
  drop_na() %>% 
  arrange(desc(abs(logFC))) %>% 
  #multiple probe IDs matched up to same gene symbol// keep one with highest abs log FC
  distinct(SYMBOL, .keep_all = TRUE) %>%  
  dplyr::slice(1:10)

write.csv(top10_fib1, 'top10_fib1.csv')

top10_fib2 <- stats_fb2 %>% 
  drop_na() %>% 
  arrange(desc(abs(logFC))) %>%  
  distinct(SYMBOL, .keep_all = TRUE) %>%  
  dplyr::slice(1:10)

write.csv(top10_fib2, 'top10_fib2.csv')

top10_epi <- stats_epi %>%  
  drop_na() %>% 
  arrange(desc(abs(logFC))) %>%  
  distinct(SYMBOL, .keep_all = TRUE) %>%  
  dplyr::slice(1:10)

write.csv(top10_epi, 'top10_epi.csv')

#create heat maps to explore the top 10 DE genes between two conditions, old and young
pdf('heatmaptop10_fib1.pdf', width = 10, height = 6)
coolmap(exprs(fibroblast1)[rownames(top10_fib1),], labRow=top10_fib1$SYMBOL, labCol = pData(fibroblast1)$condition)
dev.off()

pdf('heatmaptop10_fib2.pdf', width = 10, height = 6)
coolmap(exprs(fibroblast2)[rownames(top10_fib2),], labRow=top10_fib2$SYMBOL, labCol = pData(fibroblast2)$condition)
dev.off()

pdf('heatmaptop10_epi.pdf', width = 10, height = 6)
coolmap(exprs(epidermal)[rownames(top10_epi),], labRow=top10_epi$SYMBOL, labCol = pData(epidermal)$condition)
dev.off()


#both fib1 and fib2 are fibroblast data so let's compare the DE genes 
library(VennDiagram)
venn.diagram(x=list(top10_fib1$SYMBOL, top10_fib2$SYMBOL), category.names = c('fib1', 'fib2'),
             filename='fib1vfib2top10.png')

#top 10 DE genes in common between the two : Bcat1, Aldh1l2, Creb3l3, Ldhb
Reduce(intersect, list(top10_fib1$SYMBOL, top10_fib2$SYMBOL))



