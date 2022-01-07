\name{cDNA_synth}
\alias{cDNA_synth}
\title{
Generate cDNA synthesis cocktails
}
\description{
Generates cocktails for cDNA synthesis by reverse transcription specific to the Bio-Rad iScript cDNA Synthesis Kit. Total volume per reaction tube is 20uL.
}
\usage{
cDNA_synth(data, mass)
}
\arguments{
  \item{data}{accepts a dataframe with variables Sample and Concentration (Concentration is case-sensitive), where Sample is a character string and Concentration is a double}
  \item{mass}{accepts an numerical indicating the mass of isolated RNA (in ng)}
}
\details{
Output table is specific to the Bio-Rad iScript cDNA Synthesis Kit, and all numerical quantities are volumes in microlitres (uL).
}
\value{
  \item{Sample}{sample id carried over from initial dataframe}
  \item{Concentration}{concentration of RNA product from RNA extraction carried over from initial dataframe}
  \item{RNA_uL}{volume required to attain desired mass of isolated RNA based on Concentration}
  \item{iScript_5X_uL}{volume of iScript_5X_uL}
  \item{RTase_uL}{volume of iScript Reverse Transcriptase}
  \item{MQ_H2O_uL}{volume of Nuclease-free water or ultra-pure Milli-Q water}
  \item{total_uL}{total reaction volume}
}
\references{
https://www.bio-rad.com/sites/default/files/webroot/web/pdf/lsr/literature/4106228.pdf

https://www.emdmillipore.com/US/en/products/water-purification/type-1/5cSb.qB.89EAAAFAJqxkiQpx,nav
}
