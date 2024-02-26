read_rds_files <- function(folder_path) {
    # Get a list of all .rds files in the folder
    file_list <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)

    # Read each .rds file and save it as a global variable with the file prefix as the variable name
    for (file_path in file_list) {
        file_name <- basename(file_path)
        var_name <- gsub("\\.rds$", "", file_name)
        assign(var_name, readRDS(file_path), envir = .GlobalEnv)
    }
}
