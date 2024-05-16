#' Convert vcfR to 012 format
#'
#' @param vcf vcf file
#'
#' @return vcf file with new GT format
#' @export
#'
#' @examples
vcfR_to_012 <- function(vcf){
  library(vcfR)
  vcf <- read.vcfR(vcf)
  gt <- extract.gt(vcf, "GT")
  gt2a <- apply(gt,2, function(x) gsub("1[/|]1","2",x))
  gt2b <- gsub("0[/|]0","0",gt2a)
  gt2c <- gsub("[10][/|][10]","1",gt2b)
  gt <- gt2c

  vcf@gt <- gt
  return(vcf)
}

