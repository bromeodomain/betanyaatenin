\name{findF2}
\alias{findF2}
\title{
Find F2
}
\description{
Finds F2, the point of maximal distance from F1 (point of fusion), indicating extreme posterior (tail tip).
}
\usage{
findF2(F1, markup_fids_c1)
}
\arguments{
  \item{F1}{a vector indicating the coordinates of the point of fusion in the RAS convention}
  \item{markup_fids_c1}{a data frame including the coordinates of the fiducials defining curve C1, the exterior border of the posterior neuropore}
}
\details{
Outputs a vector defining the coordinates of F2
}
\examples{
findF2(F1, markup_fids_c1)
}
