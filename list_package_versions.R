#!/usr/bin/env Rscript
# Script to extract all R packages used in this project and their versions

# Get all R files in the project
r_files <- list.files(
  path = getwd(),
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

# Extract all library() and require() calls
all_packages <- character()

for (file in r_files) {
  lines <- readLines(file, warn = FALSE)
  # Extract library() calls
  lib_lines <- grep("^library\\(|^require\\(", lines, value = TRUE)
  for (line in lib_lines) {
    # Extract package name
    pkg <- gsub("^library\\(|^require\\(", "", line)
    pkg <- gsub("\\).*$", "", pkg)
    pkg <- gsub('"', "", pkg)
    pkg <- gsub("'", "", pkg)
    pkg <- gsub(",.*$", "", pkg)  # Remove any additional arguments like "quietly = TRUE"
    pkg <- trimws(pkg)
    if (nchar(pkg) > 0) {
      all_packages <- c(all_packages, pkg)
    }
  }
}

# Get unique packages and sort
unique_packages <- sort(unique(all_packages))

# Function to get package version
get_package_version <- function(pkg_name) {
  tryCatch({
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      pkg_info <- packageDescription(pkg_name)
      if (!is.null(pkg_info$Version)) {
        return(pkg_info$Version)
      } else {
        return("Version not found")
      }
    } else {
      return("Not installed")
    }
  }, error = function(e) {
    return("Error checking version")
  })
}

# Get versions for all packages
cat("Checking package versions...\n")
package_versions <- data.frame(
  Package = character(),
  Version = character(),
  stringsAsFactors = FALSE
)

for (pkg in unique_packages) {
  version <- get_package_version(pkg)
  package_versions <- rbind(package_versions, 
                           data.frame(Package = pkg, Version = version, 
                                     stringsAsFactors = FALSE))
  cat(sprintf("  %-30s %s\n", pkg, version))
}

# Save to file
write.table(package_versions, 
            file = "packages_with_versions.txt", 
            row.names = FALSE, 
            quote = FALSE,
            sep = "\t")

cat("\n\nPackage versions saved to: packages_with_versions.txt\n")

# Also create a formatted markdown table
cat("\n=== Creating markdown table ===\n")
markdown_table <- "| Package | Version |\n|---------|----------|\n"
for (i in 1:nrow(package_versions)) {
  markdown_table <- paste0(markdown_table, 
                          "| ", package_versions$Package[i], 
                          " | ", package_versions$Version[i], " |\n")
}

writeLines(markdown_table, "packages_versions_table.md")
cat("Markdown table saved to: packages_versions_table.md\n")

# Print summary
installed_count <- sum(package_versions$Version != "Not installed" & 
                      package_versions$Version != "Error checking version")
cat(sprintf("\nSummary: %d/%d packages are installed\n", 
            installed_count, length(unique_packages)))
