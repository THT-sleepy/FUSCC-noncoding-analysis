#导入包
library(ActiveDriverWGS)
library(data.table)

#添加fix
ADWGS_test_mod <-function (id, gr_element_coords, gr_site_coords, gr_maf, win_size, 
          this_genome, detect_depleted_mutations = FALSE) 
{
  null_res = data.frame(id, pp_element = NA, element_muts_obs = NA, 
                        element_muts_exp = NA, element_enriched = NA, pp_site = NA, 
                        site_muts_obs = NA, site_muts_exp = NA, site_enriched = NA, 
                        stringsAsFactors = F)
  gr_elements = gr_element_coords[GenomicRanges::mcols(gr_element_coords)[, 
                                                                          1] == id]
  if (length(gr_site_coords) == 0) {
    gr_sites = GenomicRanges::GRanges()
  }
  else {
    gr_sites = gr_site_coords[GenomicRanges::mcols(gr_site_coords)[, 
                                                                   1] == id]
  }
  if (length(GenomicRanges::findOverlaps(gr_elements, gr_maf)) == 
      0) {
    return(null_res)
  }
  gr_background = .create_background(gr_elements, win_size, 
                                     this_genome)
  gr_mutations = gr_maf[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_maf, 
                                                                         gr_background))]
  gr_indel = gr_mutations[GenomicRanges::mcols(gr_mutations)[, 
                                                             2] == "indel>X"]
  gr_snv = gr_mutations[GenomicRanges::mcols(gr_mutations)[, 
                                                           2] != "indel>X"]
  gr_background = GenomicRanges::setdiff(gr_background, gr_elements)
  gr_sites = GenomicRanges::intersect(gr_elements, gr_sites)
  gr_elements = GenomicRanges::setdiff(gr_elements, gr_sites)
  signt_template = .make_mut_signatures()
  trinuc_elements = .seq2signt(gr_elements, this_genome, signt_template)
  trinuc_sites = .seq2signt(gr_sites, this_genome, signt_template)
  trinuc_background = .seq2signt(gr_background, this_genome, 
                                 signt_template)
  dfr_snv = NULL
  patients_with_SNV_in_element = NULL
  if (length(gr_snv) > 0) {
    snv_elements = .mut2signt(gr_elements, gr_snv, signt_template, 
                              remove_dup_mut_per_patient = TRUE)
    snv_background = .mut2signt(gr_background, gr_snv, signt_template)
    snv_elements = merge(snv_elements, trinuc_elements, 
                         by = "signt")
    snv_background = merge(snv_background, trinuc_background, 
                           by = "signt")
    snv_elements$region = "elements"
    snv_background$region = "background"
    snv_sites = NULL
    if (length(gr_sites) > 0) {
      snv_sites = .mut2signt(gr_sites, gr_snv, signt_template, 
                             remove_dup_mut_per_patient = TRUE)
      snv_sites = merge(snv_sites, trinuc_sites, by = "signt")
      snv_sites$region = "sites"
    }
    gr_snv_el_site = gr_snv[unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_snv, 
                                                                                    c(gr_elements, gr_sites))))]
    patients_with_SNV_in_element = unique(GenomicRanges::mcols(gr_snv_el_site)[, 
                                                                               "mcols.patient"])
    dfr_snv = rbind(snv_sites, snv_elements, snv_background)
  }
  dfr_indel = NULL
  if (length(gr_indel) > 0) {
    gr_indel_fg = gr_indel[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_indel, 
                                                                            c(gr_sites, gr_elements)))]
    gr_indel_fg = gr_indel_fg[!duplicated(GenomicRanges::mcols(gr_indel_fg)[, 
                                                                            "mcols.patient"])]
    gr_indel_bg = gr_indel[S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_indel, 
                                                                            gr_background))]
    gr_indel = c(gr_indel_fg, gr_indel_bg)
    indel_tag = apply(data.frame(gr_indel)[, c("seqnames", 
                                               "start", "end", "mcols.patient")], 1, paste, collapse = "::")
    gr_indel = gr_indel[!duplicated(indel_tag)]
    indel_index_sites = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_indel, 
                                                                                gr_sites)))
    indel_index_elements = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_indel, 
                                                                                   gr_elements)))
    indel_index_background = unique(S4Vectors::queryHits(GenomicRanges::findOverlaps(gr_indel, 
                                                                                     gr_background)))
    indel_index_elements = setdiff(indel_index_elements, 
                                   indel_index_sites)
    indel_index_background = setdiff(indel_index_background, 
                                     c(indel_index_elements, indel_index_sites))
    which_indel_index_element_dup = which(GenomicRanges::mcols(gr_indel[indel_index_elements])[, 
                                                                                               "mcols.patient"] %in% patients_with_SNV_in_element)
    if (length(which_indel_index_element_dup) > 0) {
      indel_index_elements = indel_index_elements[-which_indel_index_element_dup]
    }
    which_indel_index_sites_dup = which(GenomicRanges::mcols(gr_indel[indel_index_sites])[, 
                                                                                          "mcols.patient"] %in% patients_with_SNV_in_element)
    if (length(which_indel_index_sites_dup) > 0) {
      indel_index_sites = indel_index_sites[-which_indel_index_sites_dup]
    }
    indel_sites = indel_elements = indel_background = data.frame(signt = "indel>X", 
                                                                 n_mut = NA, tri_nucleotide = "indel", n_pos = NA, 
                                                                 region = NA, stringsAsFactors = FALSE)
    indel_elements$n_mut = length(indel_index_elements)
    indel_background$n_mut = length(indel_index_background)
    indel_elements$n_pos = sum(GenomicRanges::width(gr_elements))
    indel_background$n_pos = sum(GenomicRanges::width(gr_background))
    indel_elements$region = "elements"
    indel_background$region = "background"
    if (length(gr_sites) > 0) {
      indel_sites$n_mut = length(indel_index_sites)
      indel_sites$n_pos = sum(GenomicRanges::width(gr_sites))
      indel_sites$region = "sites"
    }
    else {
      indel_sites = NULL
    }
    dfr_indel = rbind(indel_sites, indel_elements, indel_background)
  }
  dfr_mut = rbind(dfr_indel, dfr_snv)
  dfr_mut$is_site = 0 + (dfr_mut$region == "sites")
  dfr_mut$is_element = 0 + dfr_mut$region %in% c("sites", 
                                                 "elements")
  signt_with_muts = names(which(c(by(dfr_mut$n_mut, dfr_mut$signt, 
                                     sum)) > 0))
  dfr_mut = dfr_mut[dfr_mut$signt %in% signt_with_muts, , 
                    drop = FALSE]
  dfr_mut = dfr_mut[dfr_mut$n_pos > 0, , drop = FALSE]
  if(nrow(dfr_mut) > 0){
  formula_h0 = ifelse(length(signt_with_muts) > 1, "n_mut ~ signt", 
                      "n_mut ~ 1")
  h0 = stats::glm(stats::as.formula(formula_h0), offset = log(dfr_mut$n_pos), 
                  family = stats::poisson, data = dfr_mut)
  h1 = stats::update(h0, . ~ . + is_element)
  pp_element = pp_element_2way = stats::anova(h0, h1, test = "Chisq")[2, 
                                                                      5]
  coef_element = stats::coef(h1)[["is_element"]]
  element_enriched = coef_element > 0
  if (!detect_depleted_mutations & !element_enriched & !is.na(pp_element) & 
      pp_element < 0.5) {
    pp_element = 1 - pp_element_2way
  }
  if (detect_depleted_mutations & element_enriched & !is.na(pp_element) & 
      pp_element < 0.5) {
    pp_element = 1 - pp_element_2way
  }
  element_stats = .get_obs_exp(h0, dfr_mut$is_element == 1, 
                               dfr_mut, "n_mut")
  element_muts_obs = element_stats[[1]]
  element_muts_exp = element_stats[[2]]
  pp_site = site_muts_obs = site_muts_exp = site_enriched = site_depleted = NA
  if (length(gr_sites) > 0) {
    h2 = stats::update(h1, . ~ . + is_site)
    pp_site = pp_site_2way = stats::anova(h1, h2, test = "Chisq")[2, 
                                                                  5]
    coef_site = stats::coef(h2)[["is_site"]]
    site_enriched = coef_site > 0
    if (!detect_depleted_mutations & !site_enriched & !is.na(pp_site) & 
        pp_site < 0.5) {
      pp_site = 1 - pp_site_2way
    }
    if (detect_depleted_mutations & site_enriched & !is.na(pp_site) & 
        pp_site < 0.5) {
      pp_site = 1 - pp_site_2way
    }
    site_stats = .get_obs_exp(h1, dfr_mut$is_site == 1, 
                              dfr_mut, "n_mut")
    site_muts_obs = site_stats[[1]]
    site_muts_exp = site_stats[[2]]
  }
  data.frame(id, pp_element, element_muts_obs, element_muts_exp, 
             element_enriched, pp_site, site_muts_obs, site_muts_exp, 
             site_enriched, stringsAsFactors = F)}
  else{
    return(null_res)
  }
}
environment(ADWGS_test_mod) <- asNamespace('ActiveDriverWGS')
assignInNamespace("ADWGS_test", ADWGS_test_mod, ns = "ActiveDriverWGS")
args <- commandArgs(trailingOnly = TRUE)

#导入mutations
mutations <- fread(args[1])
mutations$patient <- as.character(mutations$patient)
#导入elements
all_elements_sites <- prepare_elements_from_BED4(args[2])

#run
results <- ActiveDriverWGS(mutations = mutations, elements = all_elements_sites ,ref_genome = "hg38",mc.cores=4)
#写出
write.table(results,file="results_hg38_986jiaoji.tsv",quote=F,sep="\t")

