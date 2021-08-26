library(magrittr)
library(ggplot2)

source("defs.R")

emission_prob(genotype = 0, state_tuple = c(1, 2), source_a = c(0, 2), source_d = c(2, 0))

# I have prepared a test file with the first 10,000 SNPs of the input file
dat <- readr::read_tsv("~/Data/merge_vcf.iceman_obs.wgs_short.tsv") %>%
  dplyr::transmute(
    chrom = chr,
    pos = pos,
    genotype = target_dr,
    source_a = purrr::map2(s1_anc, s2_anc, c),
    source_d = purrr::map2(s1_dr, s2_dr, c)
  )

s1_prop = 0.97
s2_prop = 1 - s1_prop

# three states in Hardy-Weinberg eq.
eq_probs <- c(s1_prop ^ 2, 2 * s1_prop * s2_prop, s2_prop ^ 2)

# transition probabilities from state i to state j are conditional on state i,
# so normalised across the second index
trans_matrix <- matrix(c(0.99, 0.01, 0.01, 0.01, 0.99, 0.01, 0.01, 0.01, 0.99), ncol=3)

fwd_dat <- run_forward(dat, eq_probs, trans_matrix, 0.1, 0.1)
fwd_vec <- fwd_dat[[1]]
scaling_facs <- fwd_dat[[2]]
bwd_vec <- run_backward(dat, scaling_facs, eq_probs, trans_matrix, 0.1, 0.1)

# posterior_normalisations:
post_norms <- colSums(fwd_vec * bwd_vec)

post_norm_matrix <- matrix(rep(post_norms, 3), nrow=3, byrow = TRUE)

posterior_probs <- fwd_vec * bwd_vec / post_norm_matrix

for_plotting <- dat %>%
  dplyr::mutate(state11 = posterior_probs[1,],
                state12 = posterior_probs[2,],
                state22 = posterior_probs[3,]) %>%
  tidyr::pivot_longer(c(state11, state12, state22), names_to = "State", values_to = "prob")

for_plotting %>%
  ggplot() + geom_area(aes(x = pos, y = prob)) + facet_grid(rows = vars(State))
