compare_svg_points <- function(old, new, tolerance){
  old <- readLines(old, warn = FALSE)
  new <- readLines(new, warn = FALSE)
  pattern <- "points='[^']*'"

  # Keep labels, styles, and SVG structure exact; tolerate only point coordinates.
  if(!identical(gsub(pattern, "points=''", old), gsub(pattern, "points=''", new)))
    return(FALSE)

  coordinates <- function(svg){
    points <- unlist(regmatches(svg, gregexpr(pattern, svg)))
    points <- gsub("points='|'", "", points)
    as.numeric(unlist(strsplit(points, "[ ,]+")))
  }
  old <- coordinates(old)
  new <- coordinates(new)
  length(old) == length(new) && all(abs(old - new) <= tolerance)
}

expect_doppelganger_tolerance <- function(title, fig, tolerance){
  testthat::local_edition(3)
  path <- tempfile(fileext = ".svg")
  on.exit(unlink(path))
  vdiffr::write_svg(fig, path, title)
  name <- paste0(gsub("^-|-$", "", gsub("[^a-z0-9]+", "-", tolower(title))), ".svg")
  testthat::expect_snapshot_file(path, name = name,
    compare = function(old, new) compare_svg_points(old, new, tolerance))
}
