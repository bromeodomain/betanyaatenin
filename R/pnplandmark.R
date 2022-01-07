reflip_coord <- function(data) {
  rl <- data[1]
  antpost <- data[2]
  sup <- data[3]

  flip_coord <- function(i) {
    if (i > 0) {
      i = i * -1
      return(i)
    } else if (i < 0) {
      i = abs(i)
      return(i)
    } else {
      return(i)
    }
  }

  R <- apply(rl, 1, flip_coord)
  A <- apply(antpost, 1, flip_coord)
  S <- sup

  data_corrected <- cbind(R, A, S)
  return(data_corrected)
}

findF2 <- function(F1, markup_fids_c1) {
  distance_f1_c1 <- function(c1) {
    axis_difference <- function(c_i, f_i) {
      if (c_i < 0 && f_i < 0) {
        diff <- abs(c_i) - abs(f_i)
      } else if (c_i > 0 && f_i > 0) {
        diff <- c_i - f_i
      } else if (c_i > 0 && f_i < 0) {
        diff <- c_i + abs(f_i)
      } else if (c_i < 0 && f_i > 0) {
        diff <- abs(c_i) + f_i
      } else {
        message("Error: invalid input")
      }
      return(diff)
    }
    rx <- axis_difference(c1[1], F1[1])
    ay <- axis_difference(c1[2], F1[2])
    sz <- axis_difference(c1[3], F1[3])
    distance = sqrt((rx^2) + (ay^2) + (sz^2))
    return(distance)
  }
  f1_c1 <- apply(markup_fids_c1, 1, distance_f1_c1)
  f1_c1_ordered <- sort(f1_c1, decreasing = TRUE)
  f2_id <- match(f1_c1_ordered[1], f1_c1)
  f2_val <- markup_fids_c1[f2_id,]
  f2 <- as.numeric(f2_val[1,])
  return(f2)
}

landmark <- function(F1, markup_fids_c1, markup_fids_c2) {
  F2 <- findF2(F1, markup_fids_c1)
  bisect1 <- rowMeans(cbind(F1, F2))
  bisect2A <- rowMeans(cbind(bisect1, F2))
  bisect2B <- rowMeans(cbind(F1, bisect1))
  v_delta <- F2-F1

  bisect_plane <- function(fids) {
    c1_on_plane <- (v_delta[1]*(fids[1]) + v_delta[2]*(fids[2]) + v_delta[3]*(fids[3]))
    return(c1_on_plane)
  }

  bisect_c1fids <- apply(markup_fids_c1, 1, bisect_plane)

  S2_constant <- (v_delta[1]*bisect1[1]) + (v_delta[2]*bisect1[2]) + (v_delta[3]*bisect1[3])
  c1_on_s2_abs <- abs(bisect_c1fids - S2_constant)
  c1_on_s2_abs_sorted <- sort(c1_on_s2_abs, decreasing = FALSE)
  c1_on_s2_closest <- match(c(c1_on_s2_abs_sorted[1], c1_on_s2_abs_sorted[2]), c1_on_s2_abs)
  F3_F4 <- markup_fids_c1[c1_on_s2_closest,]

  S6_constant <- (v_delta[1]*bisect2A[1]) + (v_delta[2]*bisect2A[2]) + (v_delta[3]*bisect2A[3])
  c1_on_s6_abs <- abs(bisect_c1fids - S6_constant)
  c1_on_s6_abs_sorted <- sort(c1_on_s6_abs, decreasing = FALSE)
  c1_on_s6_closest <- match(c(c1_on_s6_abs_sorted[1], c1_on_s6_abs_sorted[2]), c1_on_s6_abs)
  F7_F8 <- markup_fids_c1[c1_on_s6_closest,]

  S4_constant <- (v_delta[1]*bisect2B[1]) + (v_delta[2]*bisect2B[2]) + (v_delta[3]*bisect2B[3])
  bisect_c2fids <- apply(markup_fids_c2, 1, bisect_plane)
  c1_on_s4_abs <- abs(bisect_c2fids - S4_constant)
  c1_on_s4_abs_sorted <- sort(c1_on_s4_abs, decreasing = FALSE)
  c1_on_s4_closest <- match(c(c1_on_s4_abs_sorted[1], c1_on_s4_abs_sorted[2]), c1_on_s4_abs)
  F5_F6 <- markup_fids_c2[c1_on_s4_closest,]

  F1_F8 <- rbind(F1, F2, F3_F4, F5_F6, F7_F8)

  return(F1_F8)
}

secant_length <- function(F1, F2) {
  F_avg <- abs(F1 - F2)
  secant_length <- sqrt((F_avg[1])^2 + (F_avg[2])^2 + (F_avg[3])^2)
  return(secant_length)
}
