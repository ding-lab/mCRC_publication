
# set plot size in jupyter notebook
set_size <- function(w, h, factor=1.5) {
    s = 1 * factor
    options(
        repr.plot.width=w * s,
        repr.plot.height=h * s,
        repr.plot.res=100 / factor,
        jupyter.plot_mimetypes='image/png',
        jupyter.plot_scale=1
    )
}

# replace multiple mutations to "multi_hit" when plotting Oncoplot
replace_multi_hit <- function(x) {
  # Count the number of semicolons directly
  semicolon_count <- length(unlist(strsplit(x, ""))) - length(unlist(strsplit(gsub(";", "", x), "")))
  if (semicolon_count >= 2) {
    return("Multi_hit")
  } else {
    return(x)
  }
}


# remove all the  ";" 
remove_semicolons <- function(x) {
  # Use gsub to replace all semicolons with an empty string
  return(gsub(";", "", x))
}


# Function to identify and report unique values or patterns in matrix elements
report_unique_values <- function(matrix) {
  # Flatten the matrix into a vector
  elements <- as.vector(matrix) 
  # Combine all elements into a single long string, separated by semicolons
  combined_elements <- paste(elements, collapse = ";")
  # Split the combined string into individual unique values
  unique_values <- unique(unlist(strsplit(combined_elements, ";")))
  # Filter out empty strings if any
  unique_values <- unique_values[unique_values != ""]
  # Return the unique values
  return(unique_values)
}

runAllChromvar <- function(obj, assay = 'ATAC_merged') {
  DefaultAssay(obj) <- assay
  # Get a list of motif position frequency matrices from the JASPAR database
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
  )
  
  # Scan the DNA sequence of each peak for the presence of each motif
  motif.matrix <- CreateMotifMatrix(
    features = granges(obj),
    pwm = pfm,
    genome = 'BSgenome.Hsapiens.UCSC.hg38',
    use.counts = FALSE
  )
  
  # Create a new Motif object to store the results
  motif <- CreateMotifObject(
    data = motif.matrix,
    pwm = pfm
  )
  
  # Add the Motif object to the assay
  obj <- SetAssayData(
    object = obj,
    assay = assay,
    slot = 'motifs',
    new.data = motif
  )
  
  cat('doing chromvar\n')
  obj <- RegionStats(object = obj, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  obj <- RunChromVAR(
    object = obj,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  DefaultAssay(obj) <- 'chromvar'
  obj@assays$chromvar@scale.data <- obj@assays$chromvar@data
  return(obj)
}

# KellyPalette 
KellyPalette <- list(
  white = "#fdfdfd",
  black = "#1d1d1d",
  yellow = "#ebce2b",
  purple = "#702c8c",
  orange = "#db6917",
  lightBlue = "#96cde6", aqua = "#96cde6",
  red = "#ba1c30",
  buff = "#c0bd7f",
  grey = "#7f7e80", gray = "#7f7e80",
  green = "#5fa641",
  purplePink = "#d485b2", pink = "#d485b2",
  blue = "#4277b6",
  yellowPink = "#df8461", mediumOrange = "#df8461", lightBrown = "#df8461", papaya = "#df8461",
  violet = "#463397", navyBlue = "#463397", navy = "#463397",
  orangeYellow = "#e1a11a", lightOrange = "#e1a11a", manilla = "#e1a11a",
  purpleRed = "#91218c", lightPurple = "#91218c", plum = "#91218c",
  greenYellow = "#e8e948", lightYellow = "#e8e948", lemon = "#e8e948",
  redBrown = "#7e1510", brown = "#7e1510",
  yellowGreen = "#92ae31", mediumGreen = "#92ae31", lime = "#92ae31",
  yellowBrown = "#6f340d", darkBrown = "#6f340d", dirt = "#6f340d",
  redOrange = "#d32b1e", lightRed = "#d32b1e", crimson = "#d32b1e",
  oliveGreen = "#2b3514", darkGreen = "#2b3514", olive = "#2b3514"
)


theme_mydefault <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
}

