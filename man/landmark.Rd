\name{landmark}
\alias{landmark}
\title{
Landmark/outline Posterior Neuropore
}
\description{
Places 8 fiducials along curves 1 and 2, defining the posterior neuropore for geometric morphometric analysis. Begins by definining F2 using findF2 and bisecting the segment between F1 and F2 with a plane. F3 and F4 are defined as points along C1 that also satisfy the conditions of this plane. The segments created by the bisection of the segment between F1 and F2 are then themselves bisected by two more planes, defining F5 and F6 as the anterior pair along C2 and F7 and F8 as the posterior pair along C1.

Note that points must be MANUALLY REORDERED prior to executing Procrustes and PCA
}
\usage{
landmark(F1, markup_fids_c1, markup_fids_c2)
}
\arguments{
  \item{F1}{a vector indicating the coordinates of the point of fusion in the RAS convention}
  \item{markup_fids_c1}{a data frame including the coordinates of the fiducials defining curve C1, the exterior border of the posterior neuropore}
  \item{markup_fids_c2}{a data frame including the coordinates of the fiducials defining curve C2, the interior border of the posterior neuropore}
}
\details{
Outputs a table of 8 observations in the RAS convention defining the coordinates of F1 through F8
}
\examples{
landmark(F1, markup_fids_c1, markup_fids_c2)
}
