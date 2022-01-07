\name{reflip_coord}
\alias{reflip_coord}
\title{
Re-flip coordinates (LPS to RAS)
}
\description{
3DSlicer's coordinate directions (RAS) differ from the traditional DICOM coordinate system (LPS) in terms of positive and negative direction denotation. When copying and pasting values from 3DSlicer into Excel or another CSV/spreadsheet-formatting software, the L-R and P-A coordinates flip sign direction to switch from 3DSlicer's convention to the DICOM/ITK convention. reflip_coord sets these values back to the 3DSlicer convention.
}
\usage{
reflip_coord(data)
}
\arguments{
  \item{data}{accepts a data frame with three columns following the DICOM/ITK coordinate system (LPS)}
}
\details{
Outputs a three-column data frame in the 3DSlicer (RAS) coordinate system
}
\value{
  \item{R}{from left toward right}
  \item{A}{from posterior toward anterior}
  \item{S}{from inferior toward superior}
}
\references{
https://www.slicer.org/wiki/Coordinate_systems
}
\examples{
reflip_coord(data)
}
