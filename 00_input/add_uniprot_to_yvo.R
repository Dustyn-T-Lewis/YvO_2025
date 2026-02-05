#!/usr/bin/env Rscript
# Add uniprot_id and protein columns to YvO_raw.xlsx
# by mapping gene symbols to UniProt accessions

library(readxl)
library(writexl)
library(dplyr)
library(httr)
library(jsonlite)

# ---- 1. Read all three raw files ----
cat("Reading CvH_raw.xlsx...\n")
cvh <- read_excel("A_CvH_2026/00_input/CvH_raw.xlsx")
cat(sprintf("  CvH: %d rows, columns: %s\n", nrow(cvh), paste(names(cvh)[1:5], collapse=", ")))

cat("Reading HRvLR_raw.xlsx...\n")
hrlr <- read_excel("A_HRvLR_2026/00_input/HRvLR_raw.xlsx")
cat(sprintf("  HRvLR: %d rows, columns: %s\n", nrow(hrlr), paste(names(hrlr)[1:5], collapse=", ")))

cat("Reading YvO_raw.xlsx...\n")
yvo <- read_excel("A_YvO_2025/00_input/YvO_raw.xlsx")
cat(sprintf("  YvO: %d rows, columns: %s\n", nrow(yvo), paste(names(yvo)[1:3], collapse=", ")))

# ---- 2. Build gene-to-UniProt lookup from CvH and HRvLR ----
cat("\nBuilding gene -> (uniprot_id, protein) mapping from CvH and HRvLR...\n")

lookup_cvh <- cvh %>%
  select(gene, uniprot_id, protein) %>%
  distinct(gene, .keep_all = TRUE)

lookup_hrlr <- hrlr %>%
  select(gene, uniprot_id, protein) %>%
  distinct(gene, .keep_all = TRUE)

# Combine: prefer CvH first, then fill gaps with HRvLR
lookup <- bind_rows(lookup_cvh, lookup_hrlr) %>%
  distinct(gene, .keep_all = TRUE)

cat(sprintf("  Combined lookup: %d unique gene mappings\n", nrow(lookup)))

# ---- 3. Match YvO genes to lookup ----
yvo_genes <- unique(yvo$gene)
cat(sprintf("\nYvO has %d unique genes\n", length(yvo_genes)))

mapped <- yvo_genes[yvo_genes %in% lookup$gene]
unmapped <- yvo_genes[!yvo_genes %in% lookup$gene]

cat(sprintf("  Mapped from CvH/HRvLR: %d\n", length(mapped)))
cat(sprintf("  Unmapped (need UniProt query): %d\n", length(unmapped)))

# ---- 4. Query UniProt REST API for unmapped genes ----
query_uniprot_batch <- function(genes, organism = "9606") {
  results <- data.frame(
    gene = character(),
    uniprot_id = character(),
    protein = character(),
    stringsAsFactors = FALSE
  )

  # Process in batches of 25
  batch_size <- 25
  n_batches <- ceiling(length(genes) / batch_size)

  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(genes))
    batch <- genes[start_idx:end_idx]

    cat(sprintf("  Querying UniProt batch %d/%d (%d genes)...\n", i, n_batches, length(batch)))

    for (g in batch) {
      tryCatch({
        # Use UniProt REST API - search for reviewed human entries by gene name
        query_str <- sprintf('(gene_exact:"%s") AND (organism_id:%s) AND (reviewed:true)', g, organism)
        url <- sprintf(
          "https://rest.uniprot.org/uniprotkb/search?query=%s&format=json&size=1&fields=accession,gene_names,protein_name,id",
          URLencode(query_str, reserved = TRUE)
        )

        resp <- GET(url, timeout(15))

        if (status_code(resp) == 200) {
          content_json <- content(resp, as = "text", encoding = "UTF-8")
          parsed <- fromJSON(content_json, simplifyVector = FALSE)

          if (length(parsed$results) > 0) {
            entry <- parsed$results[[1]]
            acc <- entry$primaryAccession
            # protein entry name (e.g., "ALBU_HUMAN")
            entry_name <- entry$uniProtkbId
            results <- rbind(results, data.frame(
              gene = g,
              uniprot_id = acc,
              protein = entry_name,
              stringsAsFactors = FALSE
            ))
          } else {
            # Try without reviewed filter
            query_str2 <- sprintf('(gene_exact:"%s") AND (organism_id:%s)', g, organism)
            url2 <- sprintf(
              "https://rest.uniprot.org/uniprotkb/search?query=%s&format=json&size=1&fields=accession,gene_names,protein_name,id",
              URLencode(query_str2, reserved = TRUE)
            )
            resp2 <- GET(url2, timeout(15))
            if (status_code(resp2) == 200) {
              content_json2 <- content(resp2, as = "text", encoding = "UTF-8")
              parsed2 <- fromJSON(content_json2, simplifyVector = FALSE)
              if (length(parsed2$results) > 0) {
                entry2 <- parsed2$results[[1]]
                results <- rbind(results, data.frame(
                  gene = g,
                  uniprot_id = entry2$primaryAccession,
                  protein = entry2$uniProtkbId,
                  stringsAsFactors = FALSE
                ))
              } else {
                cat(sprintf("    WARNING: No UniProt match for gene '%s'\n", g))
                results <- rbind(results, data.frame(
                  gene = g,
                  uniprot_id = NA_character_,
                  protein = NA_character_,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        } else {
          cat(sprintf("    WARNING: HTTP %d for gene '%s'\n", status_code(resp), g))
          results <- rbind(results, data.frame(
            gene = g,
            uniprot_id = NA_character_,
            protein = NA_character_,
            stringsAsFactors = FALSE
          ))
        }

        Sys.sleep(0.2) # Be polite to API
      }, error = function(e) {
        cat(sprintf("    ERROR for gene '%s': %s\n", g, e$message))
        results <<- rbind(results, data.frame(
          gene = g,
          uniprot_id = NA_character_,
          protein = NA_character_,
          stringsAsFactors = FALSE
        ))
      })
    }
  }

  return(results)
}

if (length(unmapped) > 0) {
  cat("\nQuerying UniProt for unmapped genes...\n")
  uniprot_results <- query_uniprot_batch(unmapped)
  cat(sprintf("  UniProt resolved: %d / %d\n",
              sum(!is.na(uniprot_results$uniprot_id)), nrow(uniprot_results)))
} else {
  uniprot_results <- data.frame(
    gene = character(),
    uniprot_id = character(),
    protein = character(),
    stringsAsFactors = FALSE
  )
}

# ---- 5. Combine all mappings ----
full_lookup <- bind_rows(
  lookup %>% select(gene, uniprot_id, protein),
  uniprot_results
) %>%
  distinct(gene, .keep_all = TRUE)

cat(sprintf("\nFull lookup: %d entries\n", nrow(full_lookup)))

# ---- 6. Join to YvO and reorder columns ----
yvo_annotated <- yvo %>%
  left_join(full_lookup, by = "gene") %>%
  select(uniprot_id, protein, gene, description, n_seq, everything())

n_with_uniprot <- sum(!is.na(yvo_annotated$uniprot_id))
n_total <- nrow(yvo_annotated)
cat(sprintf("\nAnnotation summary:\n"))
cat(sprintf("  Total rows: %d\n", n_total))
cat(sprintf("  With UniProt ID: %d (%.1f%%)\n", n_with_uniprot, 100 * n_with_uniprot / n_total))
cat(sprintf("  Missing UniProt ID: %d\n", n_total - n_with_uniprot))

# Show missing ones if any
missing_genes <- yvo_annotated %>% filter(is.na(uniprot_id)) %>% pull(gene) %>% unique()
if (length(missing_genes) > 0) {
  cat(sprintf("  Genes without mapping: %s\n", paste(missing_genes, collapse = ", ")))
}

# ---- 7. Write output ----
output_path <- "A_YvO_2025/00_input/YvO_raw.xlsx"
cat(sprintf("\nWriting output to: %s\n", output_path))
write_xlsx(yvo_annotated, output_path)
cat("Done!\n")

# Print first few rows as verification
cat("\nFirst 5 rows of output (key columns):\n")
print(yvo_annotated %>% select(uniprot_id, protein, gene, description) %>% head(5))
