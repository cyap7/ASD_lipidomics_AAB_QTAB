#==============================================================================
#
# Metabolomics OSCA analysis
#
#==============================================================================

# Data munging prior to generating the ORMs was done in 
# `osca/metabolomics_osca_analysis.Rmd`

# Structure of this script
# 1. OREML
# - variables at the top of this file
# - make the ORMs
# - per-trait OREML + sensitivity analyses with various covariate combos
# 2. Association testing (linear models)
# - variables at the top of this file
# - per-trait OREML + sensitivity analyses with various covariate combos

data=/Volumes/YAPASD-Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/osca/data
out=/Volumes/YAPASD-Q0851/uqcyap3/ASD/Output/5_metabolomics/analysis/osca/oreml

software=~/Documents/Research/software

# Perform analyses one data type ("datan") at a time
# - used in final results of all traits EXCEPT ASD ...
datan=lipids_class_int
datan=lipids_int

# ... For the ASD analyses, where storage outliers were removed before
# calculating the ORM. This step was performed as these outliers seemed 
# to inflate the results (only for the ASD analysis). The reason for this 
# was that a subset of the ASD group were the oldest-collected samples, 
# and time since sample collection does have effects on the lipidome.
datan=lipids_class_int_rmstoragepre
datan=lipids_int_rmstoragepre

# Sensitivity analyses with different transformations of the lipidomics data
# before generating the ORMs
#datan=lipids_class_std
#datan=lipids_std
#datan=lipids_class
#datan=lipids

#------------------------------------------------------------------------------
# Variance analysis
#------------------------------------------------------------------------------

# Make the ORMs
~/osca \
--efile ${data}/${datan}.tsv \
--make-bod \
--out ${data}/${datan}

# Generate ORM + PCA
~/osca \
--befile ${data}/${datan} \
--make-orm-bin \
--orm-alg 2 \
--out ${data}/${datan}

# ~/osca \
# --befile ${data}/${datan} \
# --make-orm-gz \
# --orm-alg 2 \
# --out ${data}/${datan}

# OREML: ASD
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_asd_rmstorage_reml

# - end up using this one with "rmstoragepre" ORM (cf. post-hoc exclusion)
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemo_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemoclinical_reml

# - sensitivity analysis for above, trying to remove storage after building ORM (decided to not use this one)
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_asd_covdemo_rmstorage_reml

# - sensitivity analysis for above, including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemotime_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/asd.pheno \
#--qcovar ${data}/storage.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_asd_covstorage_reml

#---------------------------------------------
# Check macronutrient effect
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--qcovar ${data}/macronutrients.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covNdiet_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order.qcov \
--keep ${data}/macronutrients.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemoNdiet_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/asd.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_asd_covdemodiet_reml
#---------------------------------------------

# OREML: age
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/age.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemo_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemoclinical_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/age.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_age_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_age_covdemo_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemotime_reml

# OREML: Tanner
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemo_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemoclinical_reml

# - sensitivity analysis including age
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemoinclage_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/tannergenital.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_tannergenital_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_tannergenital_covdemo_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemotime_reml

# OREML: CSHQ
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemo_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemoclinical_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/sleep.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_sleep_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_sleep_covdemo_rmstorage_reml

# - sensitivity analysis including time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemotime_reml

# OREML: IQ/DQ
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemo_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemoclinical_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/iqdq.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_iqdq_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_iqdq_covdemo_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemotime_reml

# OREML: VABS motor scale
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsms.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsms_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsms.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsms_covdemo_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsms.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_vabsms_covdemo_rmstorage_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsms.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsms_covdemodiet_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsms.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsms_covdemotime_reml

# OREML: VABS gross motor
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsgm.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsgm_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsgm.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsgm_covdemo_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsgm.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsgm_covdemoclinical_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsgm.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_vabsgm_covdemo_rmstorage_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsgm.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsgm_covdemodiet_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/vabsgm.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_vabsgm_covdemotime_reml

# OREML: BMI Z-score
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_covdemo_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/bmi.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_bmi_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_bmiz_covdemo_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_covdemotime_reml

# OREML: BMI
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_covdemo_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/bmi.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_bmi_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_bmi_covdemo_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_covdemotime_reml

# OREML: sex
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch_asd.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemoasd_reml

# Error: the X^t * V^-1 * X matrix is not invertible. Please check the covariate(s) and/or the environmental factor(s).
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemoclinical_reml

# - sensitivity analysis including diet
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/age.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_age_covdemodiet_reml

# - sensitivity analysis removing storage time outliers post-ORM construction
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_sex_covdemo_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemotime_reml

# OREML: protein
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/protein.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_protein_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/protein.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_protein_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/protein.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_energy.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_protein_covdemoenergy_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/protein.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_energy.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_protein_covdemoenergy_rmstorage_reml

# OREML: fats
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/fats.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_fats_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/fats.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_fats_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/fats.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_energy.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_fats_covdemoenergy_reml

# OREML: carbohydrate
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/carbohydrate.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_carbohydrate_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/carbohydrate.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_carbohydrate_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/carbohydrate.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_energy.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_carbohydrate_covdemoenergy_reml

# OREML: sugars
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sugars.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sugars_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sugars.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sugars_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/sugars.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_energy.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sugars_covdemoenergy_reml

# OREML: cholesterol
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/cholesterol.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_cholesterol_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/cholesterol.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_cholesterol_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/cholesterol.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_energy.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_cholesterol_covdemoenergy_reml

# OREML: BSC
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bristolstool.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bristolstool_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bristolstool.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bristolstool_covdemo_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bristolstool.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bristolstool_covdemodiet_reml

#~/osca \
#--reml \
#--orm ${data}/${datan} \
#--pheno ${data}/bristolstool.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_storage_macronutrients.qcov \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_bristolstool_covdemodiet_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bristolstool.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_bristolstool_covdemodiet_rmstorage_reml

# - sensitivity analysis including collection time
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/bristolstool.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bristolstool_covdemotime_reml

# OREML: dietary PC1
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietPC1.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietPC1_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietPC1.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietPC1_covdemo_reml

# OREML: dietary PC2
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietPC2.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietPC2_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietPC2.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietPC2_covdemo_reml

# OREML: dietary PC3
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietPC3.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietPC3_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietPC3.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietPC3_covdemo_reml

# OREML: dietary diversity
~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietdiversity.pheno \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietdiversity_reml

~/osca \
--reml \
--orm ${data}/${datan} \
--pheno ${data}/dietdiversity.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--remove ${data}/outliers.id \
--out ${out}/${datan}_dietdiversity_covdemo_reml

#------------------------------------------------------------------------------
# Association analysis
#------------------------------------------------------------------------------

data=/Volumes/YAPASD-Q0851/uqcyap3/ASD/Data/5_metabolomics/analysis/osca/data
out=/Volumes/YAPASD-Q0851/uqcyap3/ASD/Output/5_metabolomics/analysis/osca/assoc

datan=lipids_class_int_rmstoragepre
datan=lipids_int_rmstoragepre
datan=lipids_class_int
datan=lipids_int
#datan=lipids_class_std
#datan=lipids_std
#datan=lipids_class
#datan=lipids

# ASD
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemo

datan=lipids_class_int_rmstoragepre
datan=lipids_int_rmstoragepre
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemoclinical

# Check effect of diet
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_dietPC.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemodiet

#~/osca \
#--befile ${data}/${datan} \
#--pheno ${data}/asd.pheno \
#--covar ${data}/sex_batch.cov \
#--qcovar ${data}/age_age2_order_storage_macronutrients.qcov \
#--logistic \
#--remove ${data}/outliers.id \
#--out ${out}/${datan}_asd_covdemodiet

# Check effect of dietary PC3 specifically
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--qcovar ${data}/dietPC3.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdietPC3only

# Check if just due to power issues - it could well be ...
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--qcovar ${data}/age_age2_order_n263dietPC3.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemo_n263dietPC3only

# Check effect of collection time
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemotime

# MOA
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_asd_rmstorage

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_asd_covdemo_rmstorage

datan=lipids_class_int_rmstoragepre
datan=lipids_int_rmstoragepre
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/asd.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_asd_covdemo

# Age
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_clinicallipids.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_time.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemotime

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch_asd.cov \
--qcovar ${data}/order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemoasd

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_age_covdemo_rmstorage

datan=lipids_class_int_rmstoragepre
datan=lipids_int_rmstoragepre
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/age.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_age_covdemo

# Sex
# ~/osca \
# --befile ${data}/${datan} \
# --pheno ${data}/sex.pheno \
# --logistic \
# --remove ${data}/outliers.id \
# --out ${out}/${datan}_sex
# 
# ~/osca \
# --befile ${data}/${datan} \
# --pheno ${data}/sex.pheno \
# --covar ${data}/batch_asd.cov \
# --qcovar ${data}/order.qcov \
# --logistic \
# --remove ${data}/outliers.id \
# --out ${out}/${datan}_sex_covdemoasd
# 
# ~/osca \
# --befile ${data}/${datan} \
# --pheno ${data}/sex.pheno \
# --orm ${data}/${datan} \
# --moa \
# --remove ${data}/outliers.id \
# --out ${out}/${datan}_sex

# Tanner
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_clinicallipids.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage_time.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemotime

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_tannergenital_covdemo_rmstorage

datan=lipids_class_int_rmstoragepre
datan=lipids_int_rmstoragepre
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/tannergenital.pheno \
--covar ${data}/sex_batch_asd.cov \
--qcovar ${data}/order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_tannergenital_covdemo

# ADOS2/G
# ~/osca \
# --befile ${data}/${datan} \
# --pheno ${data}/ados2g.pheno \
# --linear \
# --remove ${data}/outliers.id \
# --out ${out}/${datan}_ados2g
# 
# ~/osca \
# --befile ${data}/${datan} \
# --pheno ${data}/ados2g.pheno \
# --covar ${data}/sex_batch.cov \
# --qcovar ${data}/age_order.qcov \
# --linear \
# --remove ${data}/outliers.id \
# --out ${out}/${datan}_ados2g_covdemo
# 
# ~/osca \
# --befile ${data}/${datan} \
# --pheno ${data}/ados2g.pheno \
# --orm ${data}/${datan} \
# --moa \
# --remove ${data}/outliers.id \
# --out ${out}/${datan}_ados2g

# IQ/DQ
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemotime

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--orm ${data}/${datan} \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--moa \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_iqdq_covdemo_rmstorage

datan=lipids_class_int_rmstoragepre
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/iqdq.pheno \
--orm ${data}/${datan} \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_iqdq_covdemo

# CSHQ
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemotime

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers_storage.id \
--out ${out}/${datan}_sleep_covdemo_rmstorage

datan=lipids_class_int_rmstoragepre
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sleep.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--orm ${data}/${datan} \
--moa \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sleep_covdemo

# BMI z-score
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage_clinicallipids.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/order_storage_time.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmiz_covdemotime

# BMI
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmiz.pheno \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/bmi.pheno \
--covar ${data}/sex_batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--linear \
--remove ${data}/outliers.id \
--out ${out}/${datan}_bmi_covdemotime

# Sex
~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sex.pheno \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemo

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage_clinicallipids.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemoclinical

~/osca \
--befile ${data}/${datan} \
--pheno ${data}/sex.pheno \
--covar ${data}/batch.cov \
--qcovar ${data}/age_age2_order_storage_time.qcov \
--logistic \
--remove ${data}/outliers.id \
--out ${out}/${datan}_sex_covdemotime
