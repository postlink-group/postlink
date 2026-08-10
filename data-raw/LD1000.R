## code to prepare `LD1000` dataset

library(RecordLinkage)
library(fastLink)
library(dplyr)
library(clue)
library(usethis)

set.seed(2026)

### Base Data Preparation
data(RLdata10000, package = "RecordLinkage")
RLdata10000$True_ID <- identity.RLdata10000

# Format columns for fastLink
RLdata10000 <- RLdata10000 %>%
 mutate(across(c(fname_c1, lname_c1), ~as.character(ifelse(is.na(.), "", .))),
        across(c(by, bm, bd), ~as.numeric(ifelse(is.na(.), 0, .))))

# Isolate duplicate pairs
id_counts <- table(RLdata10000$True_ID)
ids_with_dups <- as.numeric(names(id_counts[id_counts == 2]))
ids_singletons <- as.numeric(names(id_counts[id_counts == 1]))

dups <- RLdata10000 %>% dplyr::filter(True_ID %in% ids_with_dups) %>%
 dplyr::arrange(True_ID)
inst1 <- dups %>% group_by(True_ID) %>% slice(1) %>% ungroup()
inst2 <- dups %>% group_by(True_ID) %>% slice(2) %>% ungroup()

# Define exact name matches vs. contains typos/mistakes
jw_f <- jarowinkler(inst1$fname_c1, inst2$fname_c1)
jw_l <- jarowinkler(inst1$lname_c1, inst2$lname_c1)
safe_flag <- (jw_f == 1 & jw_l == 1)

safe_ids <- inst1$True_ID[safe_flag]
typo_ids <- inst1$True_ID[!safe_flag]

# Sample 700 entities that can form true links
n_typos <- min(630, length(typo_ids))
n_safe <- 700 - n_typos
sampled_safe <- sample(safe_ids, n_safe)
sampled_typo <- sample(typo_ids, n_typos)
ids_true_links <- c(sampled_safe, sampled_typo)

# Add 300 unique entities to each file for a baseline mismatch rate
ids_wrong_A <- sample(ids_singletons, 300)
ids_wrong_B <- sample(setdiff(ids_singletons, ids_wrong_A), 300)

FileA <- RLdata10000 %>% dplyr::filter(True_ID %in%
                                        c(ids_true_links, ids_wrong_A)) %>%
 group_by(True_ID) %>% slice(1) %>% ungroup()

FileB_true <- RLdata10000 %>% dplyr::filter(True_ID %in% ids_true_links) %>%
 group_by(True_ID) %>% slice(2) %>% ungroup()
FileB_wrong <- RLdata10000 %>% dplyr::filter(True_ID %in% ids_wrong_B) %>%
 group_by(True_ID) %>% slice(1) %>% ungroup()

FileB <- bind_rows(FileB_true, FileB_wrong)
FileB <- FileB[sample(1:nrow(FileB)), ]

### Simulate Scientific Covariates
FileA$BMI <- scale(rnorm(1000, mean = 25, sd = 4))[,1]
FileA$Age <- scale(runif(1000, 30, 80))[,1]
FileA$Treatment <- rbinom(1000, size = 1, prob = 0.5)

# True Logistic Outcome Model
true_outcome <- data.frame(True_ID = FileA$True_ID)
eta <- -1.0 + 1.2 * FileA$BMI + 0.8 * FileA$Age - 1.5 * FileA$Treatment
true_outcome$Disease_Status <- rbinom(1000, size = 1, prob = plogis(eta))

FileB <- left_join(FileB, true_outcome, by = "True_ID")
FileB$Disease_Status[is.na(FileB$Disease_Status)] <- rbinom(300, 1,
                                                            mean(true_outcome$Disease_Status))
Oracle_A_IDs <- FileA$True_ID
Oracle_B_IDs <- FileB$True_ID

### Probabilistic Record Linkage
# Use a low threshold to generate posterior match probabilities
fl_out <- suppressMessages(fastLink(
 dfA = FileA, dfB = FileB,
 varnames = c("fname_c1", "lname_c1", "by", "bm", "bd"),
 stringdist.match = c("fname_c1", "lname_c1"),
 numeric.match = c("by", "bm", "bd"),
 threshold.match = 0.001,
 return.df = FALSE
))

### Bipartite 1:1 Matching via Linear Sum Assignment Problem (LSAP)
prob_matrix <- matrix(0, nrow = nrow(FileA), ncol = nrow(FileB))
prob_matrix[cbind(fl_out$matches$inds.a, fl_out$matches$inds.b)] <- fl_out$posterior
cost_matrix <- 1 - prob_matrix

# Solve for the optimal 1:1 assignment
assignment <- clue::solve_LSAP(cost_matrix)
Z.est <- as.numeric(assignment)

### Linkage Paradata and Latent Match Status
linked_A <- FileA[, c("BMI", "Age", "Treatment", "fname_c1", "lname_c1", "bm")] %>%
 rename(fname_A = fname_c1, lname_A = lname_c1)
linked_B <- FileB[Z.est, c("Disease_Status", "fname_c1", "lname_c1")] %>%
 rename(fname_B = fname_c1, lname_B = lname_c1)

LD1000 <- bind_cols(linked_A, linked_B)
LD1000$is_match <- Oracle_A_IDs == Oracle_B_IDs[Z.est]

# The fastLink posterior match probability for the assigned pairs
LD1000$fl_prob <- prob_matrix[cbind(1:nrow(FileA), Z.est)]

# Generate the logical safe match indicator based on exact name agreement
# and high overall RL confidence (also consider birth date)
jw_fname_post <- jarowinkler(LD1000$fname_A, LD1000$fname_B)
jw_lname_post <- jarowinkler(LD1000$lname_A, LD1000$lname_B)
LD1000$safe_match <- (jw_fname_post == 1 & jw_lname_post == 1 & LD1000$fl_prob > 0.95)

### Simulate a Clerical Review (Audit)
blocks <- as.numeric(LD1000$bm)
unique_blocks <- sort(unique(blocks))

# Randomly sample up to 20 records per birth-month block to estimate mismatch rates
m.rate_vec <- sapply(unique_blocks, function(b) {
 block_idx <- which(blocks == b)
 M_q <- length(block_idx)             # Total records in the block
 samp_size <- min(20, M_q)            # Audit sample size (m_q)

 audit_samp <- sample(block_idx, samp_size)

 # empirical correct match rate (lambda_q)
 lambda_raw <- mean(LD1000$is_match[audit_samp])

 # Based on Chambers (2009) to avoid exactly value of 1 or 0
 lambda_hat <- min((samp_size - 0.5) / samp_size, max(1 / M_q, lambda_raw))

 # Convert correct match rate to mismatch rate for the postlink package
 1 - lambda_hat
})
names(m.rate_vec) <- unique_blocks
LD1000$block_mrate_var <- m.rate_vec[as.character(LD1000$bm)]

### Linked Dataset for Secondary Analysis
# Remove direct identifiers, retain only covariates and linkage paradata
LD1000 <- LD1000 %>%
 select(Disease_Status, BMI, Age, Treatment, fl_prob, safe_match,
        bm, block_mrate_var, is_match)

# Save to package
usethis::use_data(LD1000, overwrite = TRUE)
