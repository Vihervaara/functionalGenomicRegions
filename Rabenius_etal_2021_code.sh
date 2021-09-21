#====================
#====================


wget -c -O hg38.refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

gunzip hg38.refGene.txt.gz


Rscript Rabenius_etal_2021_steps_4to6.R


bedtools intersect -v -a dREGcalls_hg38_K562.bed -b hg38_refGenes_TSSpm500.txt > enhancers.bed
bedtools intersect -u -wa -a hg38_refGenes_TSSpm500.txt -b dREGcalls_hg38_K562.bed > activeGenes_hg38_K562.bed


Rscript Rabenius_etal_2021_steps_8to12.R


echo retaining 3prime most coordinate of the bed file
awk '$6 == "+"' PROseq_K562_hg38.bed > tempPL.bed
awk '$6 == "-"' PROseq_K562_hg38.bed > tempMN.bed

awk '{$2 = $3; print}' tempPL.bed > tempPL_3p.bed
awk '{$3 = $2; print}' tempMN.bed > tempMN_3p.bed

cat tempPL_3p.bed tempMN_3p.bed | tr ' '  '\t' > temp_3p.bed
    sortBed -i temp_3p.bed > PROseq_K562_hg38_3pnt.bed

    rm *temp*



##counting engaged Pol II at promoter-proximal regions
bedtools intersect -u -wa -a PROseq_K562_hg38_3pnt.bed -b ppPolII.txt > PROseq_K562_ppPolII.bed
bedtools intersect -v -a PROseq_K562_hg38_3pnt.bed -b ppPolII.txt > ppRemoved.bed

##counting engaged Pol II at the sites of divergent transcription
bedtools intersect -u -wa -a ppRemoved.bed -b divTx.txt > PROseq_K562_ppDiv.bed
bedtools intersect -v -a ppRemoved.bed -b divTx.txt > ppdivRemoved.bed

##counting engaged Pol II at enhancers
bedtools intersect -u -wa -a ppdivRemoved.bed -b enhancers.bed > PROseq_K562_enhancers.bed
bedtools intersect -v -a ppdivRemoved.bed -b enhancers.bed > ppdivEnhRemoved.bed

##counting engaged Pol II at CPS
bedtools intersect -u -wa -a ppdivEnhRemoved.bed -b CPS.txt > PROseq_K562_CPS.bed
bedtools intersect -v -a ppdivEnhRemoved.bed -b CPS.txt > ppdivEnhCPSRemoved.bed

##counting engaged Pol II at GB
bedtools intersect -u -wa -a ppdivEnhCPSRemoved.bed -b geneBody.txt > PROseq_K562_GB.bed
bedtools intersect -v -a ppdivEnhCPSRemoved.bed -b geneBody.txt > ppdivEnhCPSgbRemoved.bed

##counting engaged Pol II at termination windows
bedtools intersect -u -wa -a ppdivEnhCPSgbRemoved.bed -b TW.txt > PROseq_K562_TW.bed
bedtools intersect -v -a ppdivEnhCPSgbRemoved.bed -b TW.txt > PROseq_K562_noGene_noEnh.bed

rm *Removed.bed



#script counts_at_functional_regions.txt


## Total count of active sites of transcription (uniquely mapping reads) in the dataset:
echo total reads:
totC=$(awk 'END { print NR }' PROseq_K562_hg38_3pnt.bed)
echo ${totC}
## Engaged Pol II at transcribed enhancers:
echo enhancers:
enhC=$(awk 'END { print NR }' PROseq_K562_enhancers.bed)
echo ${enhC}

## Engaged Pol II at divergent transcripts:
echo divergentTx:
divC=$(awk 'END { print NR }' PROseq_K562_ppDiv.bed)
echo ${divC}

## Engaged Pol II at promoters:
echo promoter:
ppC=$(awk 'END { print NR }' PROseq_K562_ppPolII.bed)
echo ${ppC}

## Engaged Pol II at gene bodies:
echo gene body:
GBC=$(awk 'END { print NR }' PROseq_K562_GB.bed)
echo ${GBC}

## Engaged Pol II at the cleavage and polyadenylation regions:
echo CPS:
CPSC=$(awk 'END { print NR }' PROseq_K562_CPS.bed)
echo ${CPSC}

## Engaged Pol II at the termination windows:
echo termination window:
twC=$(awk 'END { print NR }' PROseq_K562_TW.bed)
echo ${twC}

## Engaged Pol II at sites that did not localize to any of the identified functional regions:
echo noGenesNoEnhancers:
noGdnoE=$(awk 'END { print NR }' PROseq_K562_noGene_noEnh.bed)
echo ${noGdnoE}

#press control + D

bedtools intersect -c -wa -a ppPolII.txt -b PROseq_K562_ppPolII.bed > ppPolCounts.tmp
bedtools intersect -c -wa -a divTx.txt -b PROseq_K562_ppDiv.bed > ppDivCounts.tmp
bedtools intersect -c -wa -a enhancers.bed -b PROseq_K562_enhancers.bed > enhancerCounts.tmp
bedtools intersect -c -wa -a geneBody.txt -b PROseq_K562_GB.bed > geneBodyCounts.tmp
bedtools intersect -c -wa -a CPS.txt -b PROseq_K562_CPS.bed > CPSCounts.tmp
bedtools intersect -c -wa -a TW.txt -b PROseq_K562_TW.bed > TerminationWinCounts.tmp

awk -F '\t' -v OFS='\t' '{ $(NF+1) ="243,132,0"; print }' ppPolCounts.tmp > ppPolCounts.bed
awk -F '\t' -v OFS='\t' '{ $(NF+1) ="178,59,212"; print }' ppDivCounts.tmp > ppDivCounts.bed
awk -F '\t' -v OFS='\t' '{ $(NF+1) ="115,212,122"; print }' enhancerCounts.tmp > enhancerCounts.bed
awk -F '\t' -v OFS='\t' '{ $(NF+1) ="0,0,0"; print }' geneBodyCounts.tmp > geneBodyCounts.bed
awk -F '\t' -v OFS='\t' '{ $(NF+1) ="103,200,249"; print }' CPSCounts.tmp > CPSCounts.bed
awk -F '\t' -v OFS='\t' '{ $(NF+1) ="255,54,98"; print }' TerminationWinCounts.tmp > TerminationWinCounts.bed

cat ppPolCounts.bed ppDivCounts.bed enhancerCounts.bed geneBodyCounts.bed CPSCounts.bed TerminationWinCounts.bed > catRegions.temp

awk -F '\t' -v OFS='\t' '{ $(NF+1) ="."; print }' catRegions.temp > catRegions2.temp

awk '{print $1 "\t" $2 "\t" $3 "\t" $7 "\t" $5 "\t" $6 "\t" $2 "\t" $3 "\t" $8}' catRegions2.temp > catRegions3.temp

awk '!seen[$1,$2,$3,$6]++' catRegions3.temp | sortBed > catRegions4.temp

touch headerLine.txt
echo track name="functional_genomic_regions" itemRgb="On" >> headerLine.txt

cat headerLine.txt catRegions4.temp > functionalGenomicRegions.bed

rm *.temp
rm headerLine.txt





for x in "TBP" "GATA1" "CTCF" "H3K36me3"

do
## Factor-derived reads at promoter-proximal regions
bedtools intersect -u -wa -a K562_${x}_summits.bed -b ppPolII.txt > ${x}_K562_ppPolII.bed
bedtools intersect -v -a K562_${x}_summits.bed -b ppPolII.txt > ${x}_ppRemoved.bed

## Factor-derived reads at the sites of divergent transcription
bedtools intersect -u -wa -a ${x}_ppRemoved.bed -b divTx.txt > ${x}_K562_ppDiv.bed
bedtools intersect -v -a ${x}_ppRemoved.bed -b divTx.txt > ${x}_ppdivRemoved.bed

## Factor-derived reads at enhancers
bedtools intersect -u -wa -a ${x}_ppdivRemoved.bed -b enhancers.bed > ${x}_K562_enhancers.bed
bedtools intersect -v -a ${x}_ppdivRemoved.bed -b enhancers.bed > ${x}_ppdivEnhRemoved.bed

## Factor-derived reads at CPS
bedtools intersect -u -wa -a ${x}_ppdivEnhRemoved.bed -b CPS.txt > ${x}_K562_CPS.bed
bedtools intersect -v -a ${x}_ppdivEnhRemoved.bed -b CPS.txt > ${x}_ppdivEnhCPSRemoved.bed

## Factor-derived reads at GB
bedtools intersect -u -wa -a ${x}_ppdivEnhCPSRemoved.bed -b geneBody.txt > ${x}_K562_GB.bed
bedtools intersect -v -a ${x}_ppdivEnhCPSRemoved.bed -b geneBody.txt > ${x}_ppdivEnhCPSgbRemoved.bed

## Factor-derived reads at termination windows
bedtools intersect -u -wa -a ${x}_ppdivEnhCPSgbRemoved.bed -b TW.txt > ${x}_K562_TW.bed
bedtools intersect -v -a ${x}_ppdivEnhCPSgbRemoved.bed -b TW.txt > ${x}_K562_noGene_noEnh.bed

rm *Removed.bed
done



#script Factor_counts_at_functional_regions.txt

for x in "TBP" "GATA1" "CTCF" "H3K36me3"
do

echo total peaks of ${x}:
totC=$(awk 'END { print NR }' K562_${x}_summits.bed)
echo ${totC}

echo enhancer peaks of ${x}:
enhC=$(awk 'END { print NR }' ${x}_K562_enhancers.bed)
echo ${enhC}

echo divergentTx peaks of ${x}:
divC=$(awk 'END { print NR }' ${x}_K562_ppDiv.bed)
echo ${divC}

echo promoter peaks of ${x}:
ppC=$(awk 'END { print NR }' ${x}_K562_ppPolII.bed)
echo ${ppC}

echo genebody peaks of ${x}:
GBC=$(awk 'END { print NR }' ${x}_K562_GB.bed)
echo ${GBC}

echo CPS peaks of ${x}:
CPSC=$(awk 'END { print NR }' ${x}_K562_CPS.bed)
echo ${CPSC}

echo postCPS peaks of ${x}:
twC=$(awk 'END { print NR }' ${x}_K562_TW.bed)
echo ${twC}

echo noGenesNoEnhancers peaks of ${x}:
noGdnoE=$(awk 'END { print NR }' ${x}_K562_noGene_noEnh.bed)
echo ${noGdnoE}

done

# press control + D





##########################

