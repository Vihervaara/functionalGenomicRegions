#====================
#====================

# R script for steps 8 to 12 in Rabenius et al., 2021


Active = read.table("activeGenes_hg38_K562.bed")
refGene = read.table("hg38_refGene_allTranscripts_withHeader.txt", header=T)

refGeneAct = subset(refGene, txID %in% Active[,4])


write.table(refGeneAct[,c("chr","PPs","PPe","geneName","txID","strand")], file="ppPolII.txt", col.names=F, row.names=F, quote=F, sep="\t")

write.table(refGeneAct[,c("chr","DIVs","DIVe","geneName","txID","strand")], file="divTx.txt", col.names=F, row.names=F, quote=F, sep="\t")


shortGenes = subset(refGene, txEnd-txStart<=750)
refGeneAct_ = subset(refGeneAct, txEnd-txStart>750)



write.table(refGeneAct_[,c("chr","CPSs","CPSe","geneName","txID","strand")], file="CPS.txt", col.names=F, row.names=F, quote=F, sep="\t")

write.table(refGeneAct_[,c("chr","TWs","TWe","geneName","txID","strand")], file="TW.txt", col.names=F, row.names=F, quote=F, sep="\t")

refGeneAct_ = subset(refGeneAct_,GBe-GBs>1)  #ensuring no negative gene body lengths remain.
write.table(refGeneAct_[,c("chr","GBs","GBe","geneName","txID","strand")], file="geneBody.txt", col.names=F, row.names=F, quote=F, sep="\t")

save.image()
q()
y

