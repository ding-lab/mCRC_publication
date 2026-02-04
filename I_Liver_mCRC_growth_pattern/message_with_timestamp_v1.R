# 2024-09-18
message_with_timestamp <- function(...) {
  # Get the current timestamp
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Combine the timestamp with the original message
  msg <- paste0("[", timestamp, "] ", paste(..., collapse = " "))
  
  # Use the original message function to print the result
  base::message(msg)
}

# Set it globally by masking the base message function
assign("message", message_with_timestamp, envir = .GlobalEnv)

message("Replacing the base message function with a custom one with timestamp.")