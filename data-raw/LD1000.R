## code to prepare `LD1000` dataset

library(RecordLinkage)
library(fastLink)
library(dplyr)
library(tidyr)
library(clue)
library(usethis)

set.seed(1000)

### Base Data Preparation
data(RLdata10000, package = "RecordLinkage")
RLdata10000$True_ID <- identity.RLdata10000

# Format columns
RLdata10000 <- RLdata10000 %>%
 mutate(across(c(fname_c1, lname_c1), ~as.character(ifelse(is.na(.), "", .))),
        across(c(by, bm, bd), as.numeric)) %>%
 dplyr::filter(!is.na(bm), bm > 0, !is.na(by), by > 0)

# Isolate duplicate pairs and singles
id_counts <- table(RLdata10000$True_ID)
ids_with_dups <- as.numeric(names(id_counts[id_counts == 2]))
ids_single <- as.numeric(names(id_counts[id_counts == 1]))

dups <- RLdata10000 %>% dplyr::filter(True_ID %in% ids_with_dups) %>%
 dplyr::arrange(True_ID)
inst1 <- dups %>% group_by(True_ID) %>% slice(1) %>% ungroup()
inst2 <- dups %>% group_by(True_ID) %>% slice(2) %>% ungroup()

# Only sample true matches from the pool where birth month (bm) is identical
valid_true_ids <- inst1$True_ID[inst1$bm == inst2$bm]
inst1_valid <- inst1 %>% dplyr::filter(True_ID %in% valid_true_ids)

# Find records with relatively low degree of string typos for the true links pool
jw_f <- jarowinkler(inst1_valid$fname_c1, inst2 %>% dplyr::filter(True_ID %in% valid_true_ids) %>% pull(fname_c1))
jw_l <- jarowinkler(inst1_valid$lname_c1, inst2 %>% dplyr::filter(True_ID %in% valid_true_ids) %>% pull(lname_c1))
safe_flag <- (jw_f == 1 & jw_l == 1)

safe_ids <- inst1_valid$True_ID[safe_flag]
typo_ids <- inst1_valid$True_ID[!safe_flag]

# Sample 800 entities that can form true links
n_typos <- min(630, length(typo_ids))
n_safe <- 800 - n_typos
ids_true_links <- c(sample(safe_ids, n_safe), sample(typo_ids, n_typos))

# We take unlinked singles and group them by birth month. We then randomly
# pair two distinct people born in the same month for mismatch.
singles <- RLdata10000 %>% dplyr::filter(True_ID %in% ids_single)

singles_pairs <- singles %>%
 group_by(bm) %>%
 mutate(pair_id = rep(1:ceiling(n()/2), each = 2, length.out = n())) %>%
 group_by(bm, pair_id) %>%
 dplyr::filter(n() == 2) %>%
 ungroup()

sampled_pairs <- singles_pairs %>%
 group_by(bm, pair_id) %>%
 nest() %>%
 ungroup() %>%
 sample_n(200) %>%
 unnest(cols = c(data)) %>%
 group_by(bm, pair_id) %>%
 mutate(member = row_number()) %>%
 ungroup()

FileA_wrong <- sampled_pairs %>% dplyr::filter(member == 1)
FileB_wrong <- sampled_pairs %>% dplyr::filter(member == 2)

FileA_true <- inst1 %>% dplyr::filter(True_ID %in% ids_true_links)
FileB_true <- inst2 %>% dplyr::filter(True_ID %in% ids_true_links)

FileA <- bind_rows(FileA_true, FileA_wrong) %>% sample_frac(1) # Shuffle
FileB <- bind_rows(FileB_true, FileB_wrong) %>% sample_frac(1) # Shuffle

### Simulate Covariates
n_A <- nrow(FileA)
FileA$BMI <- round(rnorm(n_A, mean = 25, sd = 4), 1)
FileA$Age <- sample(30:80, n_A, replace = TRUE)
FileA$Treatment <- rbinom(n_A, size = 1, prob = 0.5)

# logistic outcome model
true_outcome <- data.frame(True_ID = FileA$True_ID)
eta <- -1.0 + 1.2 * scale(FileA$BMI)[,1] + 0.8 * scale(FileA$Age)[,1] - 1.5 * FileA$Treatment
true_outcome$Disease_Status <- rbinom(n_A, size = 1, prob = plogis(eta))

FileB <- left_join(FileB, true_outcome, by = "True_ID")
FileB$Disease_Status[is.na(FileB$Disease_Status)] <- rbinom(sum(is.na(FileB$Disease_Status)), 1, mean(true_outcome$Disease_Status))

Oracle_A_IDs <- FileA$True_ID
Oracle_B_IDs <- FileB$True_ID

### Probabilistic Record Linkage & 1:1 Matching
# Link on first and last name string distances.
Z.est <- numeric(nrow(FileA))
prob_matrix <- matrix(0, nrow = nrow(FileA), ncol = nrow(FileB))

unique_bms <- sort(unique(FileA$bm))

for (b in unique_bms) {
 inds_A <- which(FileA$bm == b)
 inds_B <- which(FileB$bm == b)

 if (length(inds_A) > 0 && length(inds_A) == length(inds_B)) {

  # fastLink operates on the two string comparators
  fl_out <- tryCatch({
   suppressMessages(fastLink(
    dfA = as.data.frame(FileA[inds_A, ]),
    dfB = as.data.frame(FileB[inds_B, ]),
    varnames = c("fname_c1", "lname_c1"),
    stringdist.match = c("fname_c1", "lname_c1"),
    threshold.match = 0.001,
    n.cores = 1,
    return.df = FALSE
   ))
  }, error = function(e) NULL)

  prob_b <- matrix(0, nrow = length(inds_A), ncol = length(inds_B))

  if (!is.null(fl_out) && !is.null(fl_out$matches) && nrow(fl_out$matches) > 0) {
   prob_b[cbind(fl_out$matches$inds.a, fl_out$matches$inds.b)] <- fl_out$posterior
   prob_matrix[cbind(inds_A[fl_out$matches$inds.a], inds_B[fl_out$matches$inds.b])] <- fl_out$posterior
  }

  cost_b <- 1 - prob_b
  assign_b <- clue::solve_LSAP(cost_b)
  Z.est[inds_A] <- inds_B[as.numeric(assign_b)]
 }
}

### Linkage Paradata and Latent Match Status
# Include 'by' in the merge so we can use it to define safe matches
linked_A <- FileA[, c("BMI", "Age", "Treatment", "fname_c1", "lname_c1", "bm", "by")] %>%
 rename(fname_A = fname_c1, lname_A = lname_c1, by_A = by)
linked_A$True_Disease_Status <- true_outcome$Disease_Status

linked_B <- FileB[Z.est, c("Disease_Status", "fname_c1", "lname_c1", "by")] %>%
 rename(fname_B = fname_c1, lname_B = lname_c1, by_B = by)

LD1000 <- bind_cols(linked_A, linked_B)
LD1000$is_match <- Oracle_A_IDs == Oracle_B_IDs[Z.est]
LD1000$fl_prob <- prob_matrix[cbind(1:nrow(FileA), Z.est)]

# Uses names (which were linked) and birth year
jw_fname_post <- jarowinkler(LD1000$fname_A, LD1000$fname_B)
jw_lname_post <- jarowinkler(LD1000$lname_A, LD1000$lname_B)

LD1000$safe_match <- (jw_fname_post == 1 & jw_lname_post == 1 & LD1000$by_A == LD1000$by_B)

### Simulate a Clerical Review (Audit)
blocks <- as.numeric(LD1000$bm)
unique_blocks <- sort(unique(blocks[!is.na(blocks)]))

m.rate_vec <- sapply(unique_blocks, function(b) {
 block_idx <- which(blocks == b)
 M_q <- length(block_idx)
 samp_size <- min(20, M_q)
 audit_samp <- sample(block_idx, samp_size)
 # empirical correct match rate
 lambda_raw <- mean(LD1000$is_match[audit_samp])
 # Based on Chambers (2009) to avoid exactly 0 or 1 value
 lambda_hat <- min((samp_size - 0.5) / samp_size, max(1 / M_q, lambda_raw))
 1 - lambda_hat
})
names(m.rate_vec) <- unique_blocks
LD1000$block_mrate <- m.rate_vec[as.character(LD1000$bm)]

### Linked Dataset for Secondary Analysis
LD1000 <- LD1000 %>%
 select(Disease_Status, BMI, Age, Treatment, fl_prob,
        safe_match, bm, block_mrate, True_Disease_Status, is_match)

usethis::use_data(LD1000, overwrite = TRUE)
