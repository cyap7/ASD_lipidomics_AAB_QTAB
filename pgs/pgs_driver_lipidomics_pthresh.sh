#==============================================================================
#
# PGS (P+T) for lipidomics (on UQ RCC)
#
#==============================================================================

sumstats=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/pgs/sumstats_busselton_full
pthresh=${sumstats}/pthresh
aab_qtab=/QRISdata/Q0851/uqcyap3/ASD/Data/2_GWAS/combine_aab_qtab
ref=/QRISdata/Q0851/uqcyap3/ASD/Data/2_GWAS/prs
scripts=/home/uqcyap3/scripts
#mldm=~/shared-gandalm/brain_CTP/Data/genotyping/ref/ldm/band_ukb_10k_hm3/ukb10k.mldm
pgs=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/pgs/pgs_lipids_EUR_full
snplist=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/pgs/AAB_QTAB_ABCD_sharedrsid_snplist.bim

clump=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/pgs/sumstats_clump

#==============================================================================
# Munge sumstats
#==============================================================================

#------------------------------------------------------------------------------
# 1. Intersect GWAS and target sample SNPs
#------------------------------------------------------------------------------

echo 1. BEFOREHAND: need to intersect GWAS, HM3 and target sample SNPs

#------------------------------------------------------------------------------
# 2. Clump GWAS SNPs using plink2
#------------------------------------------------------------------------------

##### TO DO: re-make list of dupsnp INCLUDING the SNP IDs eg. rs1234_A and rs1234

plink \
--bfile ${aab_qtab}/AAB_QTAB_MERGED_1_22_ACRC_sex_ref_info0.3_maf01_hwe1e6 \
--make-bed \
--exclude ${aab_qtab}/AAB_QTAB_MERGED_1_22_ACRC_sex_ref_info0.3_maf01_hwe1e6.snpdup \
--out ${aab_qtab}/AAB_QTAB_MERGED_1_22_ACRC_sex_ref_info0.3_maf01_hwe1e6_rmmulti

for filen in PC.O.18.0.20.4. PE.P.19.0.20.4...b.
do
plink \
--bfile ${aab_qtab}/AAB_QTAB_MERGED_1_22_ACRC_sex_ref_info0.3_maf01_hwe1e6_rmmulti \
--clump ${sumstats}/${filen}.sumstats_AAB_QTAB_ABCD_rmdup \
--keep ${clump}/lipids_id.txt \
--clump-snp-field SNP \
--clump-field P \
--clump-p1 0.5 \
--clump-p2 0.5 \
--clump-r2 0.10 \
--clump-kb 250 \
--out ${pthresh}/${filen}_AAB_QTAB_ABCD_r2clump0.1.sumstats
done

#------------------------------------------------------------------------------
# 3. Intersect GWAS and clumped SNPs
#------------------------------------------------------------------------------

for filen in PC.O.18.0.20.4. PE.P.19.0.20.4...b.
do
qsub -v scripts=${scripts},\
f1=${sumstats}/${filen}.sumstats_AAB_QTAB_ABCD_rmdup,\
c1=1,\
h1="header1",\
f2=${pthresh}/${filen}_AAB_QTAB_ABCD_r2clump0.1.sumstats.clumped,\
c2=3,\
h2="header2",\
mc="SNP",\
ot="separate",\
od=${pthresh},\
n1=${filen}_AAB_QTAB_ABCD_r2clump0.1.sumstats \
${scripts}/intersect_2files.sh
done

#------------------------------------------------------------------------------
# 4. Threshold ${gwas} at various cut-offs
#------------------------------------------------------------------------------

# https://choishingwan.github.io/PRS-Tutorial/plink/
for filen in PC.O.18.0.20.4. PE.P.19.0.20.4...b.
do
# for --q-score-range: generate snp.pvalue
awk '{print $1,$7}' ${pthresh}/${filen}_AAB_QTAB_ABCD_r2clump0.1.sumstats > ${pthresh}/snp.pvalue_${filen}

# for --q-score-range: generate range_list
echo "5e-08 0 5e-08" > ${pthresh}/range_list_${filen}
echo "5e-07 0 5e-07" >> ${pthresh}/range_list_${filen}
echo "5e-06 0 5e-06" >> ${pthresh}/range_list_${filen}
echo "5e-05 0 5e-05" >> ${pthresh}/range_list_${filen}
echo "5e-04 0 5e-04" >> ${pthresh}/range_list_${filen}
echo "0.005 0 0.005" >> ${pthresh}/range_list_${filen}
done

#------------------------------------------------------------------------------
# 5. Generate PGS at each of these cut-offs
#------------------------------------------------------------------------------

target=AAB_QTAB
for filen in PC.O.18.0.20.4. PE.P.19.0.20.4...b.
do
~/plink2 \
--bfile ${aab_qtab}/AAB_QTAB_MERGED_1_22_ACRC_sex_ref_info0.3_maf01_hwe1e6 \
--score ${pthresh}/${filen}_AAB_QTAB_ABCD_r2clump0.1.sumstats 1 2 5 list-variants \
--q-score-range ${pthresh}/range_list_${filen} ${pthresh}/snp.pvalue_${filen} \
--out ${pthresh}/${target}_pgs.${filen}.pthresh
done
