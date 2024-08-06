# -------------------------------------------------------------------------------------------------- #
# Description: Data download, pre-processing and select TNBC samples
# -------------------------------------------------------------------------------------------------- #
# Load packages and create directory
# -------------------------------------------------------------------------------------------------- #

setwd(".")
dir.base <- getwd()
source(file.path(dir.base, "code/utility/utility_function.R"))
dir.exp <- file.path(dir.base, "/Data/Expr/")
dir.exp.counts <- file.path(dir.base, "/Data/Expr/counts")
dir.cData <- file.path(dir.base, "/Data/cData/")
dir.pcr <- file.path(dir.cData, "pCR")
dir.raw <- file.path(dir.base, "/Data/raw/")
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)


# -------------------------------------------------------------------------------------------------- #
# GSE25055
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE25055", baseDir = dir.raw) 
untar(file.path(dir.raw, "GSE25055/GSE25055_RAW.tar"), exdir = file.path(dir.raw, "GSE25055/"))
affyob <- ReadAffy(celfile.path = file.path(dir.raw, "GSE25055/"))
pinfo_25055_65 <- pData(phenoData(getGEO("GSE25066")[[1]]))
TNBC_sample <- read.csv(file.path(dir.cData, "GSE25055_GSE25065_pinfo.csv"), row.names = 1)

pinfo_25055_65c <- pinfo_25055_65[match(rownames(TNBC_sample), rownames(pinfo_25055_65)),c(56:80)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IQR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133a')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133a.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get pCR

pCR <- TNBC_sample$pcr_rd
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)]
pCR <- data.frame(pCR = pCR[names])
identical(rownames(pCR), colnames(exp2)[-1]) # TRUE
pinfo_25055c <- pinfo_25055_65c[match(colnames(exp2)[-1], rownames(pinfo_25055_65c)),]

pCR <- na.omit(pCR)

write.csv(pinfo_25055c, file.path(dir.cData, "pinfo_GSE25055c.csv"))
write.csv(pCR, file.path(dir.pcr,"/pCR_GSE25055.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE25055c.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE25065
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE25065", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE25065/GSE25065_RAW.tar"), exdir = file.path(dir.raw, "GSE25065/"))
affyob<- ReadAffy(celfile.path = file.path(dir.raw,"GSE25065/"))
pinfo_25055_65 <- pData(phenoData(getGEO("GSE25066")[[1]]))
TNBC_sample <- read.csv(file.path(dir.cData,"/GSE25055_GSE25065_pinfo.csv"), row.names = 1)
pinfo_25055_65c <- pinfo_25055_65[match(rownames(TNBC_sample), rownames(pinfo_25055_65)),c(56:80)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133a')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133a.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- TNBC_sample$pcr_rd
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)]
pCR <- data.frame(pCR = pCR[names])
identical(rownames(pCR), colnames(exp2)[-1]) # TRUE
pinfo_25065c <- pinfo_25055_65c[match(colnames(exp2)[-1], rownames(pinfo_25055_65c)),]

pCR <- na.omit(pCR)

write.csv(pinfo_25065c, file.path(dir.cData, "pinfo_GSE25065c.csv"))
write.csv(pCR, file.path(dir.pcr, "pCR_GSE25065.csv"))
write_csv(exp2, file.path(dir.exp, 'exp_GSE25065c.csv'))


# -------------------------------------------------------------------------------------------------- #
# GSE41998
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE41998", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE41998/GSE41998_RAW.tar"), exdir = file.path(dir.raw, "GSE41998/"))
affyob <- ReadAffy(celfile.path = file.path(dir.raw, "GSE41998/"))
pinfo_41998 <- pData(phenoData(getGEO("GSE41998")[[1]]))
TNBC_sample <- read.csv(file.path(dir.cData, "/GSE25066_GSE41998_allsample_meta.csv"),row.names = 1)
rownames(TNBC_sample) <- TNBC_sample$SAMPLE_ID

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp<-exprs(object) 

colnames(expp)<-substr(colnames(expp),1,10)

#### Filter TNBC samples

TNBC_sample = TNBC_sample %>% filter(TNBC == "YES")
exp_TNBC = expp[,colnames(expp) %in% TNBC_sample$SAMPLE_ID]

#### use probe with largest IDR to represent the gene

arrayIQR<-apply(exp_TNBC,1,IQR)
probe<-rownames(exp_TNBC)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,'hgu133a')

exp2<-exp_TNBC[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu133a.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR = TNBC_sample %>% mutate(pCR = ifelse(PCR == "YES", "pCR", ifelse(PCR == "NO", "RD", NA)))
pCR = pCR$pCR
names(pCR) = TNBC_sample$SAMPLE_ID
names = names(pCR)[names(pCR) %in% colnames(exp2)]
pCR = data.frame(pCR = pCR[names])
pinfo_41998c <- TNBC_sample[match(colnames(exp2)[-1], rownames(TNBC_sample)),]
identical(rownames(pinfo_41998c), colnames(exp2)[-1]) # TRUE

pCR = na.omit(pCR)

write.csv(pinfo_41998c, file.path(dir.cData, "/pinfo_GSE41998c.csv"))
write.csv(pCR, file.path(dir.pcr,"/pCR_GSE41998.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE41998c.csv'))


# -------------------------------------------------------------------------------------------------- #
# GSE16446
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE16446", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE16446/GSE16446_RAW.tar"), exdir = file.path(dir.raw, "GSE16446/"))
affyob <- ReadAffy(celfile.path = file.path(dir.raw, "GSE16446/"))
pinfo_16446 <- pData((getGEO("GSE16446")[[1]]))
pinfo_16446c <- pinfo_16446[,c(49:64)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples

array <- expp
ER <- c("205225_at")
PR <- c("208305_at")
HER2 <- c("216836_s_at")
all <- c(ER, PR, HER2)

array <- array[all,]

PR <- array[PR,]
truehist(PR,h=0.1)
lines(density(PR))
neg <- PR[PR<mean(PR)]
pos <- PR[PR>mean(PR)]
p0 <- c(p=length(neg)/length(PR), u1=mean(neg), s1=sd(neg),u2=mean(pos), s2=sd(pos))
es <- optim(p0,mix.obj, x=PR)
post2 <- es$par[1]*dnorm((PR-es$par[2])/es$par[3])/(es$par[1]*dnorm((PR-es$par[2])/es$par[3])+(1-es$par[1])*(dnorm((PR-es$par[4])/es$par[5])))

aa2 <- post2
for (i in 1:length(aa2)){
  if ( aa2[i] > 0.5 ) {aa2[i] <- 1} else
  {aa2[i] <- 0}
}

TN <- cbind(aa2)
colnames(TN) <- c("PR")
pinfo <- cbind(pinfo_16446c, TN)

TNBC_sample <- pinfo %>% filter(PR == 1, `esr1bimod:ch1` == 0,
                            `erbb2bimod:ch1`  == 0)

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133plus2')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133plus2.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- TNBC_sample %>% mutate(pCR = ifelse(`pcr:ch1` == 1, "pCR", ifelse(`pcr:ch1` == 0, "RD", NA)))
pCR <- pCR$pCR
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR <- data.frame(pCR = pCR[names])
identical(rownames(pCR), colnames(exp2)[-1]) # TRUE

pCR = na.omit(pCR)

write.csv(pinfo_16446c, file.path(dir.cData,"/pinfo_GSE16446c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE16446.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE16446c.csv'))


# -------------------------------------------------------------------------------------------------- #
# GSE18728
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE18728", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE18728/GSE18728_RAW.tar"), exdir = file.path(dir.raw, "GSE18728/"))
affyob <- ReadAffy(celfile.path = file.path(dir.raw, "GSE18728/"))
pinfo_18728 <- pData(phenoData(getGEO("GSE18728")[[1]]))
pinfo_18728c <- pinfo_18728[,c(39:47)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp<-exprs(object) 

colnames(expp)<-substr(colnames(expp),1,9)

#### Filter TNBC samples

TNBC_sample = pinfo_18728c %>% filter(`er original:ch1`  == "neg", `pr:ch1` == "neg",
                               `her2 summary:ch1` == "neg")

exp_TNBC = expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR<-apply(exp_TNBC,1,IQR)
probe<-rownames(exp_TNBC)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,'hgu133plus2')

exp2<-exp_TNBC[uniqueGenes,]
geneSymbol<-getSYMBOL(rownames(exp2),"hgu133plus2.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR = TNBC_sample %>% mutate(pCR = ifelse(`response category:ch1` == "R", "pCR", 
                                          ifelse(`response category:ch1` == "NR", "RD", NA)))
pCR = pCR$pCR
names(pCR) = rownames(TNBC_sample)
names = names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR = data.frame(pCR = pCR[names])
identical(rownames(pCR), colnames(exp2)[-1]) # TRUE

pCR = na.omit(pCR)

write.csv(pinfo_18728c,file.path(dir.cData, "/pinfo_GSE18728c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE18728.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE18728c.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE20194
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE20194", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE20194/GSE20194_RAW.tar"), exdir = file.path(dir.raw, "GSE20194/"))
affyob <- ReadAffy(celfile.path = file.path(dir.raw, "GSE20194/"))
pinfo_20194 <- pData(phenoData(getGEO("GSE20194")[[1]]))
pinfo_20194c <- pinfo_20194[,c(52:67)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples


TNBC_sample <- pinfo_20194c %>% filter(`er_status:ch1` == "N", `her2 status:ch1` == "N",
                               `pr_status:ch1` == "N")

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe<-rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133a')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133a.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- TNBC_sample$`pcr_vs_rd:ch1`
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR <- data.frame(pCR = pCR[names])

pCR <- na.omit(pCR)

write.csv(pinfo_20194c, file.path(dir.cData, "/pinfo_GSE20194c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE20194.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE20194c.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE20271
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE20271", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE20271/GSE20271_RAW.tar"), exdir = file.path(dir.raw, "GSE20271/"))
affyob <- ReadAffy(celfile.path=file.path(dir.raw, "GSE20271/"))
pinfo_20271 <- pData(phenoData(getGEO("GSE20271")[[1]]))
pinfo_20271c <- pinfo_20271[,c(58:82)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples


TNBC_sample <- pinfo_20271c %>% filter(`er status:ch1` == "N", `her 2 status:ch1` == "N",
                                       `pr status:ch1`   == "N")

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes<-findLargest(as.vector(probe),arrayIQR,'hgu133a')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133a.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- TNBC_sample$`pcr or rd:ch1`
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR <- data.frame(pCR = pCR[names])

pCR = na.omit(pCR)

write.csv(pinfo_20271c, file.path(dir.cData, "/pinfo_GSE20271c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE20271.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE20271c.csv'))


# -------------------------------------------------------------------------------------------------- #
# GSE18864
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE18864", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE18864/GSE18864_RAW.tar"), exdir = file.path(dir.raw, "GSE18864/"))
affyob<- ReadAffy(celfile.path = file.path(dir.raw, "GSE18864/"))
pinfo_18864 <- pData(phenoData(getGEO("GSE18864")[[1]]))
pinfo_18864c <- pinfo_18864[,c(37:42)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples


TNBC_sample <- pinfo_18864c %>% filter(`er/pr/her2 status:ch1` == "neg/neg/neg" & `brca genotype:ch1` != 'Unknown')

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133plus2')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133plus2.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- ifelse(TNBC_sample$`miller-payne response:ch1` %in% c("0", "1", "2", "3"), "RD", 
             ifelse(TNBC_sample$`miller-payne response:ch1` %in% c("4", "5"), "pCR", NA))
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR <- data.frame(pCR = pCR[names])

pCR <- na.omit(pCR)

write.csv(pinfo_18864c, file.path(dir.cData, "/pinfo_GSE18864c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE18864.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE18864c.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE22093
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE22093", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE22093/GSE22093_RAW.tar"), exdir = file.path(dir.raw, "GSE22093/"))
affyob<- ReadAffy(celfile.path = file.path(dir.raw, "GSE22093/"))
pinfo_22093 <- pData(phenoData(getGEO("GSE22093")[[1]]))
pinfo_22093c <- pinfo_22093[,c(50:62)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples

array <- expp
ER <- c("205225_at")
PR <- c("208305_at")
HER2 <- c("216836_s_at")
all <- c(ER, PR, HER2)

array<-array[all,]

PR <- array[PR,]
truehist(PR,h=0.1)
lines(density(PR))
neg <- PR[PR<mean(PR)]
pos <- PR[PR>mean(PR)]
p0 <- c(p=length(neg)/length(PR), u1=mean(neg), s1=sd(neg),u2=mean(pos), s2=sd(pos))
es <- optim(p0,mix.obj, x=PR)
post2 <- es$par[1]*dnorm((PR-es$par[2])/es$par[3])/(es$par[1]*dnorm((PR-es$par[2])/es$par[3])+(1-es$par[1])*(dnorm((PR-es$par[4])/es$par[5])))

HER2 <- array[HER2,]
neg <- HER2[HER2<11.5]
pos <- HER2[HER2>11.5]
p0 <- c(p=length(neg)/length(HER2), u1=mean(neg), s1=sd(neg),u2=mean(pos), s2=sd(pos))
es <- optim(p0,mix.obj, x=HER2)
truehist(HER2,h = 0.1)
lines(density(HER2))
post3 <- es$par[1]*dnorm((HER2-es$par[2])/es$par[3])/(es$par[1]*dnorm((HER2-es$par[2])/es$par[3])+(1-es$par[1])*(dnorm((HER2-es$par[4])/es$par[5])))


aa2 <- post2
for (i in 1:length(aa2)){
  if (aa2[i]>0.5 ) {aa2[i]<-1} else
  {aa2[i]<-0}
}

aa3 <- post3
for (i in 1:length(aa3)){
  if (aa3[i]>0.5 ) {aa3[i]<-1} else
  {aa3[i]<-0}
}

TN <- cbind(aa2,aa3)
colnames(TN) <- c("PR","HER2")
pinfo <- cbind(pinfo_22093c, TN)

TNBC_sample <- pinfo %>% filter(PR == 1, HER2 == 1,
                              `er positive vs negative by esr1 mrna gene expression (probe 205225_at):ch1`  == "ERneg")

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133a')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133a.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- TNBC_sample$`pcr.v.rd:ch1`
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR <- data.frame(pCR = pCR[names])

pCR = na.omit(pCR)

write.csv(pinfo_22093c, file.path(dir.cData, "/pinfo_GSE22093c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE22093.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE22093c.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE32646
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE32646", baseDir = dir.raw)
untar(file.path(dir.raw, "GSE32646/GSE32646_RAW.tar"), exdir = file.path(dir.raw, "GSE32646/"))
affyob <- ReadAffy(celfile.path = file.path(dir.raw, "GSE32646/"))
pinfo_32646 <- pData(phenoData(getGEO("GSE32646")[[1]]))
pinfo_32646c <- pinfo_32646[,c(43:52)]

#### fRMA

object <- frma(affyob, verbose=TRUE)
expp <- exprs(object) 

colnames(expp) <- substr(colnames(expp),1,9)

#### Filter TNBC samples


TNBC_sample <- pinfo_32646c %>% filter(`er status ihc:ch1` == "negative", `her2 status fish:ch1` == "negative",
                                     `pr status ihc:ch1` == "negative")

exp_TNBC <- expp[,colnames(expp) %in% rownames(TNBC_sample)]

#### use probe with largest IDR to represent the gene

arrayIQR <- apply(exp_TNBC,1,IQR)
probe <- rownames(exp_TNBC)
uniqueGenes <- findLargest(as.vector(probe),arrayIQR,'hgu133plus2')

exp2 <- exp_TNBC[uniqueGenes,]
geneSymbol <- getSYMBOL(rownames(exp2),"hgu133plus2.db")
exp2 <- exp2 %>% as.data.frame() %>% mutate(
  Symbol = geneSymbol,
  .before = 1
)

#### Get clinical info

pCR <- ifelse(TNBC_sample$`pathologic response pcr ncr:ch1` == "nCR", "RD", 
             ifelse(TNBC_sample$`pathologic response pcr ncr:ch1` == "pCR", "pCR", NA))
names(pCR) <- rownames(TNBC_sample)
names <- names(pCR)[names(pCR) %in% colnames(exp2)[-1]]
pCR <- data.frame(pCR = pCR[names])

pCR <- na.omit(pCR)

write.csv(pinfo_32646c, file.path(dir.cData, "/pinfo_GSE32646c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE32646.csv"))
write_csv(exp2, file.path(dir.exp, '/exp_GSE32646c.csv'))

# -------------------------------------------------------------------------------------------------- #
# VUMC
# -------------------------------------------------------------------------------------------------- #

#### Re-annotation RNA-seq count data

VUMC.count <- read_csv(file.path(dir.raw, "VUMC/20150630_BRE0904_tnbc_gene.count2.csv"))
mapped_ens <- readRDS(file.path(dir.raw, "VUMC/mapped_ens.rds"))
pinfo_VUMC <- read_csv(
  file.path(dir.cData, "Vanderbilt_meta.csv")
) %>% filter( Timepoint == "PRE")

exp <- VUMC.count %>% filter(Feature_gene_biotype == "protein_coding")
anno_gene <- exp[,c(1:4)]
map <- mapped_ens[mapped_ens$ensembl_gene_id %in% exp$Feature,]
exp$Feature_gene_name[which(exp$Feature %in% mapped_ens$ensembl_gene_id)] <- map$hgnc_symbol

exp2 <- exp[,colnames(exp) %in% pinfo_VUMC$RNA_seq_ID]
exp2 <- exp2 %>% mutate(
  Symbol = exp$Feature_gene_name,
  .before = 1
)
identical(colnames(exp2)[-1], pinfo_VUMC$RNA_seq_ID)
dim(exp2) 
# [1] 19898    45

write.csv(exp2, file.path(dir.exp.counts, "/VUMC_counts.csv"))
rm(exp)

#### Re-annotation RNA-seq fpkm data

VUMC_fpkm <- read_csv(file.path(dir.raw, "VUMC/20150630_BRE0904_tnbc_gene.fpkm.rmdup.csv"))

gene_keep <- intersect(anno_gene$Feature_gene_name, VUMC_fpkm$Feature_gene_name)
exp <- VUMC_fpkm[VUMC_fpkm$Feature_gene_name %in% gene_keep,]
exp <- cbind(anno_gene, exp[match(anno_gene$Feature_gene_name, exp$Feature_gene_name),])
exp <- exp[,-2]
map <- mapped_ens[mapped_ens$ensembl_gene_id %in% exp$Feature,]
exp$Feature_gene_name[which(exp$Feature_gene_name %in% mapped_ens$ensembl_gene_id)] = map$hgnc_symbol
exp2 <- exp[, colnames(exp) %in% pinfo_VUMC$RNA_seq_ID]
exp2 <- exp2 %>% mutate(
  Symbol = exp$Feature_gene_name,
  .before = 1
)

VUMC.fpkm2 <- exp2[rowSums(exp2[,-1]) > 0,]
min <- min(VUMC.fpkm2[,-1][VUMC.fpkm2[,-1]>0])
VUMC.fpkm2[,-1][VUMC.fpkm2[,-1]<=0] <- VUMC.fpkm2[,-1][VUMC.fpkm2[,-1]<=0]+min

exp2 <- data.frame(Symbol = VUMC.fpkm2$Symbol, log2(VUMC.fpkm2[,-1]))

pCR <- pinfo_VUMC %>% mutate(pCR = ifelse(pCR == "NO", "RD", "pCR"), .keep = "none")

write.csv(pCR, file.path(dir.pcr, "/pCR_VUMC.csv"))
write_csv(exp2, file.path(dir.exp, "/VUMC_fpkm_log2.csv"))

# -------------------------------------------------------------------------------------------------- #
# GSE22226
# -------------------------------------------------------------------------------------------------- #

#### Get data

gset <- getGEO("GSE22226", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
fvarLabels(gset) <- make.names(fvarLabels(gset))

exp <- exprs(gset)
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  exp <- log2(ex) 
}
exp <- data.frame(Symbol = gset@featureData@data[["Gene.symbol"]], exp, stringsAsFactors = F) %>%
  filter(Symbol != "") %>% na.omit() 
exp <- exp[!duplicated(exp$Symbol),]

pinfo_22226 <- pData(phenoData(gset))[,c(2, 60:80)]

rownames(pinfo_22226) <- pinfo_22226$geo_accession

#### TNBC sample

TNBC_sample <- pinfo_22226 %>% filter(`er (0=negative; 1=positive;):ch2` == 0, 
                                     `her2 (0=negative; 1=positive;):ch2` == 0, 
                                     `pgr  (0=negative; 1=positive;):ch2` == 0)
clin_info <- data.frame(
  Age = TNBC_sample$`age:ch2`, 
  pCR = TNBC_sample$`pathological complete response (pcr):ch2`,
  Grade = TNBC_sample$`histologic grade (1=grade i (low);  2= grade ii (intermediate); 3= grade iii (high); 4=indeterminate):ch2`,
  event = TNBC_sample$`relapse-free survival indicator (1=event; local or distant progression or death, 0=censor at last follow-up):ch2`,
  time = as.numeric(gsub("[^[:digit:].]", "", TNBC_sample$`relapse-free survival time – time from chemo start date until earliest:ch2`))/365)

clin_info <- clin_info %>% mutate(pCR = ifelse(pCR == "Yes", 1, ifelse(pCR == "No", 0, NA)))
rownames(clin_info) <- rownames(TNBC_sample)
clin_info <- na.omit(clin_info)

exp <- exp[,c(1,match( rownames(clin_info), colnames(exp)))]

write.csv(pinfo_22226, file.path(dir.cData, "/pinfo_GSE22226c.csv"))
write.csv(clin_info, file.path(dir.pcr, "/clin_GSE22226.csv"))
write_csv(exp, file.path(dir.exp,"/exp_GSE22226.csv"))

# -------------------------------------------------------------------------------------------------- #
# GSE149322
# -------------------------------------------------------------------------------------------------- #

#### Get data

gset = getGEO("GSE149322", GSEMatrix =TRUE, AnnotGPL=TRUE)

exp <- lapply(gset, exprs)

pinfo_149322 <- do.call(rbind, lapply(
    gset, function(g) pData(phenoData(g))
  )
)

pinfo_149322c <- pinfo_149322[,c(36:39)]
rownames(pinfo_149322c) <- pinfo_149322$geo_accession

exp <- lapply(
  gset, function(g){
    e <- data.frame(Symbol = g@featureData@data[,3], exprs(g)) %>%
      filter(Symbol != "") %>% na.omit() 
    #colnames(e) <- e[,1]
    e
  }
) 

intersect_gene <- lapply(exp, function(e) e$Symbol) %>% purrr::reduce(intersect)
exp2 <- lapply(exp, function(e) e[match(intersect_gene, e$Symbol),]) %>% purrr::reduce(merge, by = "Symbol")


#### TNBC sample

TNBC_sample <- pinfo_149322c %>% filter(`her2:ch1` == 0, `hr:ch1` == 0)
exp_TNBC <- exp2[,c(1, match(rownames(TNBC_sample), colnames(exp2)))]

clin_info = data.frame(pCR = TNBC_sample$`pcr:ch1`, Arm = TNBC_sample$`arm:ch1`)

write.csv(pinfo_149322c, file.path(dir.cData, "/pinfo_GSE149322c.csv"))
write.csv(clin_info, file.path(dir.pcr, "/clin_info_GSE149322.csv"))
write_csv(exp_TNBC, file.path(dir.exp, "/exp_GSE149322.csv"))

# -------------------------------------------------------------------------------------------------- #
# GSE154524
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE154524", baseDir = dir.raw)
exp_154524 = read.table(
  file.path(dir.raw, "GSE154524/GSE154524_CALGB40603.389PreTreat.JUNCTION.salmon_gene_normalized.matrix_FOR_GEO.v2.txt.gz"),
  header = T
) 
colnames(exp_154524)[1] <- "Symbol"
pinfo_154524 <- read.csv(
  file.path(dir.cData, "/GSE154524_pinfo.csv"),
  row.names = 1
)
gset <- getGEO("GSE154524")
pinfo_154524_supp <- do.call(rbind, lapply(
    gset, function(g) pData(phenoData(g))
  )
)

pinfo_154524_supp$title <- gsub("\\-", "\\.", pinfo_154524_supp$title)
rownames(pinfo_154524_supp) <- pinfo_154524_supp$geo_accession
pinfo_154524_supp <- pinfo_154524_supp[match(pinfo_154524$PID, pinfo_154524_supp$title), ]
rownames(pinfo_154524) <- pinfo_154524_supp$geo_accession

exp <- exp_154524[,c(1,match(pinfo_154524$PID, colnames(exp_154524)))]
colnames(exp)[-1] <- rownames(pinfo_154524)

#### log2 transform

exp_154524_log2 <- data.frame(Symbol = exp$Symbol,log2(exp[,-1] + 1))

pCR <- pinfo_154524 %>% mutate(pCR = ifelse(pCR_Breast == 1, "pCR", ifelse(pCR_Breast == 0, "RD", NA)), .keep = "none")
identical(rownames(pCR), colnames(exp)[-1])
pCR <- na.omit(pCR)

write.csv(pinfo_154524, file.path(dir.cData, "/pinfo_GSE154524c.csv"))
write.csv(pCR, file.path(dir.pcr, "/pCR_GSE154524.csv"))
write_csv(exp_154524_log2, file.path(dir.exp, '/exp_GSE154524_log2.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE192341
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE192341",  baseDir = dir.raw)
gset <- getGEO("GSE192341", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

exp <- data.table::fread(file.path(dir.raw, "GSE192341/GSE192341_processed_data.txt.gz"), header = T) %>%
  column_to_rownames("ensembl_gene_id")
exp <- exp[exp$V93 != "",]
info <- exp[,c(88:93)]
exp2 <- exp[,-c(88:93)]
exp2 <- exp2 %>% dplyr::mutate(
  Symbol = info$V93,
  .before = 1
)

pinfo_192341 <- pData(phenoData(gset))
pinfo_192341c <- pinfo_192341[,c(2, 44:53)]

exp2 <- exp2[,c(1,match(pinfo_192341c$`sample_id:ch1`, colnames(exp2)))]
colnames(exp2)[-1] <- pinfo_192341c$geo_accession

#### TNBC 

TNBC_sample <- pinfo_192341c %>% filter(`subtype:ch1` %in% 'TN')
exp2_TN <- exp2[,colnames(exp2) %in% c("Symbol", TNBC_sample$geo_accession)]
identical(colnames(exp2_TN)[-1], TNBC_sample$geo_accession)
pCR <- TNBC_sample %>% mutate(pCR = ifelse(`pathological complete response (pcr):ch1` == 'pCR', "pCR", 
                      ifelse(`pathological complete response (pcr):ch1` == 'No pCR', "RD",NA)) ,.keep = "none")

write.csv(pinfo_192341c,  file.path(dir.cData, "/pinfo_GSE192341c.csv"))
write_csv(exp2_TN, file.path(dir.exp.counts, '/exp_GSE192341_counts.csv'))
write.csv(pCR,  file.path(dir.pcr, "//pCR_GSE192341.csv"))

#### count to tpm

f_length <- (info$V92 - info$V91)
tpms <- apply(exp2_TN[,-1], 2, function(x) round(tpm(x, f_length), 4))

#### log2 transform

tpms_log2 <- data.frame(Symbol = exp2_TN$Symbol, log2(tpms + 1))

write_csv(tpms_log2, file.path(dir.exp,'/exp_GSE192341_tpms_log2.csv'))

# -------------------------------------------------------------------------------------------------- #
# GSE163882
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE163882",  baseDir = dir.raw)
gset <- getGEO("GSE163882", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]
exp_163882 <- read.csv(file.path(dir.raw, "GSE163882/GSE163882_all.data.tpms_222Samples.csv.gz"), row.names = 1)
pinfo_163882 <- pData(phenoData(gset))

pinfo_163882c <- pinfo_163882[,c(23, 48:55)]

TN_163882 = pinfo_163882c %>% filter(`estrogen receptor status:ch1` == 'N' & 
                                     `her2 receptor status:ch1` == 'N' & 
                                     `progesterone receptor status:ch1` == 'N')


TN_exp <- exp_163882[, match(TN_163882$description, colnames(exp_163882))]
colnames(TN_exp) <- rownames(TN_163882)

TN_exp_log2 <- log2(TN_exp + 1)
TN_exp_log2 <- TN_exp_log2 %>% mutate(
  Symbol = exp_163882$annotation,
  .before = 1
)

pCR <- TN_163882 %>% mutate(pCR = `response to nac:ch1`, .keep = "none")
identical(rownames(pCR), colnames(TN_exp_log2)[-1])

write.csv(pinfo_163882c,  file.path(dir.cData, "/pinfo_GSE163882c.csv"))
write_csv(TN_exp_log2, file.path(dir.exp, '/exp_GSE163882_tpms_log2.csv'))
write.csv(pCR, file.path(dir.pcr, "pCR_GSE163882.csv"))

# -------------------------------------------------------------------------------------------------- #
# GSE164458
# -------------------------------------------------------------------------------------------------- #

#### Get data

getGEOSuppFiles("GSE164458",  baseDir = dir.raw)
exp_164458 <- read.table(
  file.path(dir.raw, "GSE164458/GSE164458_BrighTNess_RNAseq_log2_Processed_ASTOR.txt.gz"),
  header = T
) 
colnames(exp_164458) <- gsub("\\X", "", colnames(exp_164458))
pinfo_164458 <- read.csv(
  file.path(dir.cData, "/GSE164458_meta.csv"),
  row.names = 1
)
gset <- getGEO("GSE164458")[[1]]
pinfo_164458_supp <- pData(phenoData(gset))
pinfo_164458_supp$title <- do.call(rbind,stringr::str_split(pinfo_164458_supp$title,"_"))[,1]
pinfo_164458c <- pinfo_164458[match(rownames(pinfo_164458_supp), rownames(pinfo_164458)),]

exp_164458c <- exp_164458[,-1][,match(colnames(exp_164458)[-1], pinfo_164458_supp$title)]
colnames(exp_164458c) <- rownames(pinfo_164458_supp)
identical(colnames(exp_164458c), rownames(pinfo_164458c))
exp_164458c <- exp_164458c %>% mutate(
  Symbol = exp_164458$Sample,
  .before = 1
)

write.csv(pinfo_164458c,  file.path(dir.cData, "/pinfo_GSE164458c.csv"))
write_csv(exp_164458c, file.path(dir.exp, '/exp_GSE164458c.csv'))

