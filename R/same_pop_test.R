#' @title A simple falsification test for the Two-Sample Mendelian randomisation 'same population' assumption.
#' @description This function is a falsification test for the 'same population' assumption made by Two-Sample Mendelian randomisation. Please note that this function tests the absolute difference in SNP effects such that negative effects imply that the exposure GWAS has smaller betas than the outcome GWAS
#' @param Bexp SNP-phenotype associations for the exposure sample.
#' @param Bout SNP-phenotype associations for the outcome sample.
#' @param SEexp Standard error of the SNP-phenotype associations for the exposure sample.
#' @param SEout Standard error of the SNP-phenotype associations for the outcome sample.
#' @param SNPlist list of rsid's for the SNPs used.
#' @param Random Logical indicating if a random effects meta-analysis should run, defult = T.
#' @param Fixed Logical indicating if a fiexed effects meta-analysis should run, defult = T.
#' @param sm A character string indicating underlying summary measure, e.g., "RD", "RR", "OR", "ASD", "HR", "MD", "SMD", or "ROM".
#' @keywords 2SMR
#' @param Fisher Logical indicating if to Fisher's methods (meta-analysis of p-vlalues) to combine differences
#' @export
#' @examples
#' library(TwoSampleMR)
#' #extracting list of SNPs used as instruments
#' exp <- extract_instruments(outcomes="ieu-b-4760")
#' #extracting data on outcome from outcome sample
#' out<- extract_outcome_data( exp$SNP,  c("ieu-a-1009"))
#' #extracting data on outcome from exposure sample
#' ukb_out<-extract_outcome_data(dat$SNP,"ukb-b-4062")
#' #combining data for analysis
#' ukb_out$beta.ukb<-ukb_out$beta.outcome
#' ukb_out$se.ukb<-ukb_out$se.outcome
#' ukb_out<-ukb_out[,c("SNP","beta.ukb","se.ukb","effect_allele.outcome")]
#' snp_comp<-merge(dat[,c("SNP","beta.outcome","effect_allele.exposure","se.outcome")],ukb_out, by="SNP")
#' #running analysis
#' same_pop_test(SNPlist=snp_comp$SNP, Bexp=snp_comp$beta.ukb, Bout=snp_comp$beta.outcome, SEexp=snp_comp$se.ukb, SEout=snp_comp$se.outcome, Random=T, Fixed=F)
#' #example taken from (insert DIO)

#Imports: "meta"  ## taken from here:  https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html  n.b. this does not work htough
same_pop_test <- function(Bexp, Bout, SEexp, SEout, SNPlist, Fisher=F, Random=T, Fixed=T,sm="md"){
  Bout<-Bout*sign(Bexp)
  Bexp<-abs(Bexp)
  diff<-Bexp-Bout
  sediff<-(SEexp^2+SEout^2)^0.5
  if(Fisher==T){
    p.diff<-2*pnorm(abs((diff)/sediff),low=F)
    return(RecordTest::fisher.method(p.diff))
  }
  if(Fisher==F){
    meta::forest(meta.comp <- meta::metagen(TE=diff, seTE=sediff ,   test.subgroup=F, studlab = SNPlist, fixed = Fixed, random=Random,  sm=sm, tau.common = FALSE),  leftcols = c("studlab","effect.ci") ,leftlabs=c("SNP","Difference [95% CI]"), smlab="", rightcols =F)
    return(meta.comp) }}

