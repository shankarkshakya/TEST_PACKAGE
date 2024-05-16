# Function to process VCF file and generate CSV with unique SNPs
#' @title Find unique genotype for each sample from vcf file
#'
#' @param vcf_file a vcf file
#'
#' @return a csv file for each sample with unique genotype
#' @export
#'
#' @examples

find_unique_GT_per_sample <- function(vcf_file) {
  library(vcfR)
  library(ggplot2)
  library(forcats)

  vcf <- read.vcfR(vcf_file)
  gt <- extract.gt(vcf, "GT", return.alleles = TRUE)

  vcf <- vcf[which(!duplicated(gt)), ]
  gt <- extract.gt(vcf, "GT")

  gt2a <- apply(gt,2, function(x) gsub("1[/|]1","1",x))
  gt2b <- gsub("0[/|]0","0",gt2a)
  gt2c <- gsub("[10][/|][10]","0.5",gt2b)
  gt <- gt2c


  count_list <- apply(gt, 1, table)
  my_min_list <- vector("list", length = length(count_list))

  for(i in 1:length(count_list)) {
    ith_item <- count_list[[i]]
    my_min_list[[i]] <- min(ith_item)
  }

  q <- gt[which(unlist(lapply(my_min_list, function(x) x == 1))), ]
  q <- data.frame(q, stringsAsFactors = FALSE)

  test <- q
  my_unique <- vector("list", length = nrow(test))
  names(my_unique) <- rownames(test)


  # Extract unique SNPs per sample
  for (i in 1:nrow(test)){
    p <- data.frame(apply(test[i,], 1, table))
    index <- which(p[,1] == 1)
    my_unique[[i]] <- colnames(test[which(test[i,] == rownames(p)[index])])
  }
  # Process each series to find unique SNPs
  atleast_one_unique_gt <- test
  #rownames(atleast_one_unique_gt) <- atleast_one_unique_gt$X
  #atleast_one_unique_gt <- atleast_one_unique_gt[, -1]
  dim(atleast_one_unique_gt)

  for (i in 1:ncol(atleast_one_unique_gt)){
    ith_series <- colnames(atleast_one_unique_gt)[i]
    print(ith_series)

    keep <- names(which(unlist(my_unique) == ith_series))
    p <- atleast_one_unique_gt[rownames(atleast_one_unique_gt) %in% keep, ]

    filename <- paste0("Unique_for_", ith_series, ".csv")

    write.csv(p, filename, row.names = TRUE)
  }
}

