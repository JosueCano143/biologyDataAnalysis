> load("TCGA_COADREAD_comp_data.RData")
> View(tcga_coadread)
> view(tcga_coadread_class)
Error in view(tcga_coadread_class) : 
  no se pudo encontrar la función "view"
> View(tcga_coadread_class)
> y_index <- which(tcga_coadread_class == "Young")
> o_index <- which(tcga_coadread_class == "Old")
> young <- tcga_coadread[,y_index]
> old <- tcga_coadread[,o_index]
> View(young)
> View(old)
> p_values <- matrix(data = NA, nrow(tcga_coadread),ncol = 4)
> rownames(p_values) <- c(rownames(tcga_coadread))
> colnames(p_values) <- c("Young", "Old", "p-value", "difMean")
> View(p_values)
> for (i in c(1:nrow(p_values))) {test <- t.test(young[i,],old[i,])
+ p_values[i,] <- c(test$estimate[1], test$estimate[2], test$p.value, (test$estimate[1]-test$estimate[2]))
+ }
> order(p_values[,3])[1:10]
[1] 11267 14077  6724 10949 17702 13076  9629 13449  8813 20006
> p_values[order(p_values[,3])[1:10],]
Young       Old      p-value     difMean
MTERF        8.39193346 7.6258149 8.054799e-07  0.76611852
PRND         4.16365896 2.9448635 4.748019e-06  1.21879543
FZD9         4.51446092 3.3978873 1.000739e-05  1.11657367
MLF1         8.38633153 7.4555142 1.630553e-05  0.93081730
TBC1D3P2     0.01050615 0.0882879 6.441823e-05 -0.07778175
PCSK1N       5.16845587 3.5492721 7.661858e-05  1.61918372
LOC100009676 7.22734164 7.0113538 9.199702e-05  0.21598788
PIWIL1       4.27838391 5.7870158 9.945316e-05 -1.50863187
KCNRG        4.13745416 4.6079456 1.726251e-04 -0.47049146
ZNF239       7.58108540 8.0794047 1.856189e-04 -0.49831925
> index <- which(p_values[,1:2] < 1)
> new_p_values <- p_values[-index,]
> View(new_p_values)
> sorted_p_values <- new_p_values[order(new_p_values[,3]),]
> View(sorted_p_values)
> write.table(rownames(sorted_p_values[which(sorted_p_values[,4] > 0),])[1:10], sep="\t", quote=F, row.names=F, col.names=F)
MTERF
PRND
FZD9
MLF1
PCSK1N
LOC100009676
KCNS3
C5orf49
CCDC90B
TRMT12
> write.table(rownames(sorted_p_values[which(sorted_p_values[,4] < 0),])[1:10], sep="\t", quote=F, row.names=F, col.names=F)
PIWIL1
KCNRG
ZNF239
ZNF600
YOD1
ATG16L1
AVPI1
ANKMY1
UBE2MP1
C5orf56
> new_p_values[order(p_values[,3])[1:10],]
Error in new_p_values[order(p_values[, 3])[1:10], ] : 
  subíndice fuera de  los límites
> new_p_values[order(new_p_values[,3])[1:10],]
Young      Old      p-value    difMean
MTERF        8.391933 7.625815 8.054799e-07  0.7661185
PRND         4.163659 2.944864 4.748019e-06  1.2187954
FZD9         4.514461 3.397887 1.000739e-05  1.1165737
MLF1         8.386332 7.455514 1.630553e-05  0.9308173
PCSK1N       5.168456 3.549272 7.661858e-05  1.6191837
LOC100009676 7.227342 7.011354 9.199702e-05  0.2159879
PIWIL1       4.278384 5.787016 9.945316e-05 -1.5086319
KCNRG        4.137454 4.607946 1.726251e-04 -0.4704915
ZNF239       7.581085 8.079405 1.856189e-04 -0.4983192
KCNS3        8.485033 7.886361 1.893069e-04  0.5986718
> write.table(rownames(sorted_p_values[which(sorted_p_values[,4] > 0),])[1:50], sep="\t", quote=F, row.names=F, col.names=F)
MTERF
PRND
FZD9
MLF1
PCSK1N
LOC100009676
KCNS3
C5orf49
CCDC90B
TRMT12
GAL
SPATA17
GATA4
ZNF75A
GAMT
CCDC67
DKK1
BHMT
MYL3
TFAP2E
TPPP3
SH3GL2
HLTF
NAT8L
ACOT13
HAVCR1
FRAS1
CDH7
HSCB
DSC3
TMSB15A
RPL39L
GLT8D1
FOLR1
ZMAT5
SLC46A1
GHDC
GREB1L
LOC348926
SLC8A2
TBCK
C16orf71
LEFTY2
AHSG
MAP9
GLDC
CBY1
HSPC157
MGA
GDF11
> write.table(rownames(sorted_p_values[which(sorted_p_values[,4] < 0),])[1:50], sep="\t", quote=F, row.names=F, col.names=F)
PIWIL1
KCNRG
ZNF239
ZNF600
YOD1
ATG16L1
AVPI1
ANKMY1
UBE2MP1
C5orf56
TRIM26
DUSP5
STK39
CELSR3
PPP1R15B
SOX9
IFNG
ELF3
AGAP7
BATF2
FZD5
ANXA11
LIF
IL15RA
SLC25A28
HOXB8
ZFP36
EDAR
IRF1
TAP2
GPD2
ATF3
SOX1
ZC3H11A
YBX2
LPGAT1
NAPB
HSPA1B
FOXD4L1
CNTD2
NR4A2
ZNF384
CLK2P
SF3B4
MET
LRIT2
ZBTB7B
SOCS1
TOB1
MAT1A
> index_order_mean <- order(abs(p_values[,4]), decreasing=T)
> order_mean <- p_values[index_order_mean,]
> View(order_mean)
> index_order_p_value <- order(p_values[,3])
> order_p_value <- p_values[index_order_p_value,]
> View(order_p_value)
> rownames(order_mean)[1:20]
[1] "GATA4"        "PCSK1N"       "PIWIL1"       "XIST"         "DUSP27"       "HAVCR1"       "DSC3"        
[8] "DKK1"         "PRND"         "FOLR1"        "CPS1"         "GAL"          "SOX1"         "FZD9"        
[15] "GLDC"         "LOC100190940" "GREB1L"       "SULT1E1"      "FEZF1"        "BHMT"        
> rownames(order_p_value)[1:20]
[1] "MTERF"        "PRND"         "FZD9"         "MLF1"         "TBC1D3P2"     "PCSK1N"       "LOC100009676"
[8] "PIWIL1"       "KCNRG"        "ZNF239"       "KCNS3"        "C5orf49"      "CCDC90B"      "TRMT12"      
[15] "GAL"          "SPATA17"      "GATA4"        "ZNF75A"       "GAMT"         "CCDC67"      
> par(mfrow=c(1,2))
> plot(index_order_mean[1:100],index_order_mean[1:100], xlab = "Medias", ylab = "Medias")
> abline(a=1,b=1)
> plot(index_order_mean[1:100],index_order_p_value[1:100], xlab = "Medias", ylab = "Medias")
> abline(a=1, b=1)