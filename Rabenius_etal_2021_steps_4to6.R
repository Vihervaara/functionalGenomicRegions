#====================
#====================

# R script for steps 4 to 6 in Rabenius et al., 2021



refGene = read.table("hg38.refGene.txt", header=F)

names(refGene) =c("bin", "txID", "chr", "strand", "txStart", "txEnd","cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "geneName", "cdsStartStat", "cdsEnsStat", "exonFrames")
refGene = refGene[-grep("_", refGene$chr),]    #maintains chromosomes 1-22, X, Y and M.
refGene$chr = factor(refGene$chr) #drops the extra levels removed above.

refGene = refGene[,c(3,5:6,4,13,2)]
head(refGene) #look at the dataframe again

refGene_pl = subset(refGene, strand=="+")
refGene_mn = subset(refGene, strand=="-")


### Genes on the plus strand:

refGene_pl$TSS = refGene_pl$txStart    # TSS
refGene_pl$CPS = refGene_pl$txEnd        # CPS

refGene_pl$DIVs = refGene_pl$txStart-750    # region of divergent transcription
refGene_pl$DIVe = refGene_pl$txStart-251

refGene_pl$PPs = refGene_pl$TSS-250    # promoter-proximal region
refGene_pl$PPe = refGene_pl$TSS+249

refGene_pl$GBs = refGene_pl$TSS+250    # genebody
refGene_pl$GBe = refGene_pl$CPS-501

refGene_pl$CPSs = refGene_pl$CPS-500    # CPS region
refGene_pl$CPSe = refGene_pl$CPS+499

refGene_pl$TWs = refGene_pl$CPS+500    # termination window
refGene_pl$TWe = refGene_pl$CPS+10499


#### Genes on the minus strand:

refGene_mn$TSS = refGene_mn$txEnd        # TSS
refGene_mn$CPS = refGene_mn$txStart    # CPS

refGene_mn$DIVs = refGene_mn$txEnd+251    # divergent transcription region
refGene_mn$DIVe = refGene_mn$txEnd+750

refGene_mn$PPs = refGene_mn$TSS-249    # promoter-proximal region
refGene_mn$PPe = refGene_mn$TSS+250

refGene_mn$GBs = refGene_mn$CPS+501    # genebody
refGene_mn$GBe = refGene_mn$TSS-250

refGene_mn$CPSs = refGene_mn$CPS-499    # CPS region
refGene_mn$CPSe = refGene_mn$CPS+500

refGene_mn$TWs = refGene_mn$CPS-10499    # termination window
refGene_mn$TWe = refGene_mn$CPS-500


#### combine the data of plus and minus stands:

refGene = rbind(refGene_pl, refGene_mn)

refGene$promC1 = refGene$TSS-500
refGene$promC2 = refGene$TSS+500


#### generate data files:

write.table(refGene, file="hg38_refGene_allTranscripts.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(refGene, file="hg38_refGene_allTranscripts_withHeader.txt", col.names=T, row.names=F, quote=F, sep="\t")

write.table(refGene[,c("chr","promC1","promC2","txID","geneName")], file="hg38_refGenes_TSSpm500.txt", col.names=F, row.names=F, quote=F, sep="\t")


save.image()
q()
y    


