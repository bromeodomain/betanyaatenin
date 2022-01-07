\name{generate_qPCR}
\alias{generate_qPCR}
\title{
Generate qPCR cocktails
}
\description{
Generates cocktails for each target gene for qPCR given the number of samples, number of trials per sample, number of distinct genotypes, and the names of the target genes. Total number of samples is increased by one to account for pipetting error before calculations are performed. Total volume per reaction tube is 10uL.
}
\usage{
generate_qPCR(num_sample, trials, genotypes, genes)
}
\arguments{
  \item{num_sample}{integer value indicating the number of samples if single genotype, number of pairs if paired WT/cKO, or number of triplets if tripled WT/cKO/Het}
  \item{trials}{number of trials/tubes per num_sample}
  \item{genotypes}{number of distinct genotypes (quantity may be any integer)}
  \item{genes}{a character vector delineating the target genes to be tested}
}
\details{
Output table is specific to qPCR performed with Takara TB Green Advantage qPCR, and all numerical values in the tables are in microlitres (uL).
}
\value{
  \item{tb_green}{volume of Takara TB Green Advantage qPCR 2X premix. Final concentration is 1X}
  \item{primer}{assumes a 1:10 dilution from 100x storage stock to 10x diluted stock, and the final concentration given by generate_qPCR provides a 1X working concentration}
  \item{mq_h2o}{volume of ultra-pure Milli-Q water}
  \item{template}{voume of cDNA product from reverse transcription cDNA synthesis}
}
\references{
https://www.takarabio.com/products/real-time-pcr/real-time-pcr-kits/qpcr-with-tb-green-detection/tb-green-advantage-qpcr-premixes

https://www.emdmillipore.com/US/en/products/water-purification/type-1/5cSb.qB.89EAAAFAJqxkiQpx,nav
}
\seealso{
cDNA_synth
}
