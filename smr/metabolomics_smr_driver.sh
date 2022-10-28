#==============================================================================
#
# Lipidomics SMR 
#
#==============================================================================

sumstats_lipid_raw=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/pgs/sumstats_busselton
sumstats_lipid=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/smr/sumstats_busselton
sumstats_pheno=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/smr/sumstats_pheno
besd=/QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/smr/besd
ldref=/QRISdata/Q0851/uqcyap3/ASD/Data/2_GWAS/combine_aab_qtab/AAB_QTAB_MERGED_1_22_ACRC_sex_ref_info0.3_maf01_hwe1e6_EUR
out=/QRISdata/Q0851/uqcyap3/ASD/Output/5_metabolomics/analysis/smr

#==============================================================================
# Format lipid sumstats files --> save as .esd
#==============================================================================

#------------------------------------------------------------------------------
# Make .esd files
#------------------------------------------------------------------------------

# Lipid classes
for filen in Total.PC.O. Total.PE.O. Total.PE.P. Total.LPE.P.
do
cat <(echo "Chr SNP Bp A1 A2 Freq Beta se p") <(awk '{print $3,$2,$4,$6,$7,$8,$9,$10,$11}' ${sumstats_lipid_raw}/${filen}.txt | tail -n+2) > ${sumstats_lipid}/${filen}.esd
done

# Lipid species
for filen in PC.O.18.0.20.4. PC.O.42.4...b. PC.P.35.2...a. PC.P.35.2...b. LPC.O.22.0. LPC.P.17.0...a. PC.15.0.22.6. PC.17.0.22.6. PC.O.16.0.22.6. PC.P.18.0.22.6. PE.P.17.0.20.4...b. PE.P.17.0.22.6...a. PE.P.17.0.22.6...b. PE.P.18.0.22.6. PE.P.19.0.20.4...b. PE.P.20.0.22.6. PE.P.20.1.20.4. PE.P.20.1.22.6. dimethyl.CE.22.6. 
do
cat <(echo "Chr SNP Bp A1 A2 Freq Beta se p") <(awk '{print $4,$3,$5,$7,$8,$9,$10,$11,$12}' ${sumstats_lipid_raw}/${filen}.txt | tail -n+2) > ${sumstats_lipid}/${filen}.esd
done

#------------------------------------------------------------------------------
# Make .besd files
#------------------------------------------------------------------------------

# Manually generated .flist file
# - only included lipid classes/species that had GWAS sumstats

~/software/smr_Linux \
--eqtl-flist ${besd}/lipids_sumstats_sig_asd_iq_sleep.flist \
--make-besd \
--out ${besd}/lipids_sumstats_sig_asd_iq_sleep

# IQ
~/software/smr_Linux \
--eqtl-flist ${besd}/iq_PC.O.18.0.20.4..flist \
--make-besd \
--out ${besd}/iq_PC.O.18.0.20.4.

# Sleep
~/software/smr_Linux \
--eqtl-flist ${besd}/sleep_PE.P.19.0.20.4...b..flist \
--make-besd \
--out ${besd}/sleep_PE.P.19.0.20.4...b.

#==============================================================================
# (pheno sumstats are already formatted)
#==============================================================================

cp /QRISdata/Q0851/uqcyap3/ASD/Data/2_GWAS/prs/IQ.sumstats ${sumstats_pheno}
cp /QRISdata/Q0851/uqcyap3/ASD/Data/2_GWAS/prs/SleepDuration.sumstats ${sumstats_pheno}

#==============================================================================
# Run SMR
#==============================================================================

phenon=IQ
~/software/smr_Linux \
--bfile ${ldref} \
--gwas-summary ${sumstats_pheno}/${phenon}.sumstats \
--beqtl-summary ${besd}/iq_PC.O.18.0.20.4. \
--out ${out}/${phenon}_PC.O.18.0.20.4.

phenon=SleepDuration
~/software/smr_Linux \
--bfile ${ldref} \
--gwas-summary ${sumstats_pheno}/${phenon}.sumstats \
--beqtl-summary ${besd}/sleep_PE.P.19.0.20.4...b. \
--out ${out}/${phenon}_PE.P.19.0.20.4...b.

# Runs SMR + plots
phenon=IQ
~/software/smr_Linux \
--bfile ${ldref} \
--gwas-summary ${sumstats_pheno}/${phenon}.sumstats \
--beqtl-summary ${besd}/iq_PC.O.18.0.20.4. \
--plot \
--gene-list /QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/smr/glist-hg19 \
--probe PC.O.18.0.20.4. \
--probe-wind 500 \
--out ${out}/${phenon}_PC.O.18.0.20.4._plot

phenon=SleepDuration
~/software/smr_Linux \
--bfile ${ldref} \
--gwas-summary ${sumstats_pheno}/${phenon}.sumstats \
--beqtl-summary ${besd}/sleep_PE.P.19.0.20.4...b. \
--plot \
--gene-list /QRISdata/Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/smr/glist-hg19 \
--probe PE.P.19.0.20.4...b. \
--probe-wind 500 \
--out ${out}/plot/${phenon}_PE.P.19.0.20.4...b._plot

#==============================================================================
# LD clump
#==============================================================================

module load plink

filen=PC.O.18.0.20.4.
plink \
--bfile ${ldref} \
--clump ${sumstats_lipid}/${filen}.esd \
--clump-r2 0.5 \
--clump-snp-field SNP \
--clump-field p \
--out ${out}/${filen}_r2clump0.5.sumstats

filen=PE.P.19.0.20.4...b.
plink \
--bfile ${ldref} \
--clump ${sumstats_lipid}/${filen}.esd \
--clump-r2 0.5 \
--clump-snp-field SNP \
--clump-field p \
--out ${out}/${filen}_r2clump0.5.sumstats

