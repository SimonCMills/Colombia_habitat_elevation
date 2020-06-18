## Copy new/updated files to/from cluster 
##  note: need to run this going from local to cluster before running cluster code
##      need to run this going from cluster to local to retrieve new results from 
##      cluster
##
##  source() from script that specifies from & target variables
if(!any(ls() == "target")) stop("need to specify target")
if(!any(ls() == "from")) stop("need to specify from")

# if target directory doesn't already exist, create it
if(!file.exists(target)) dir.create(target)

# get all files in respective directories
from_files <- list.files(from, full.names = T, recursive = T)
target_files <- list.files(target, full.names = T, recursive = T)

# note fnames_toMatch has the target path
from_finfo <- file.info(from_files) %>%
    mutate(fnames_from = row.names(.), 
           fnames_toMatch = gsub(from, target, fnames_from))

target_finfo <- file.info(target_files) %>%
    mutate(fnames_toMatch = row.names(.))

# join to identify additions, removals, and changes
matched <- full_join(from_finfo, target_finfo, by="fnames_toMatch", 
                     suffix=c("_from", "_target"))

toAdd <- matched %>% 
    filter(is.na(mtime_target))

toUpdate <- matched %>% 
    filter(!is.na(mtime_target), 
           mtime_from > mtime_target)

toRemove <- matched %>%
    filter(is.na(mtime_from))

## iterate over the respective categories
# create subdirectories that don't already exist
gsub("(^.*)\\/.*$", "\\1", toAdd$fnames_toMatch) %>% unique %>%
    lapply(., dir.create, showWarnings = F)

if(nrow(toAdd) != 0) {
    apply(toAdd, 1, 
          function(x) file.copy(x["fnames_from"], 
                                to = gsub("(^.*)\\/.*$", "\\1", gsub(from, target, x["fnames_toMatch"])), 
                                recursive = T))   
}
if(nrow(toUpdate) != 0) {
    apply(toUpdate, 1, 
          function(x) file.copy(x["fnames_from"], 
                                to = gsub("(^.*)\\/.*$", "\\1", gsub(from, target, x["fnames_toMatch"])), 
                                recursive = T))
}
if(nrow(toRemove) != 0) {
    apply(toRemove, 1, 
          function(x) file.remove(x["fnames_toMatch"]))
}

print(paste0("Added ", nrow(toAdd), " file(s)"))
print(paste0("Updated ", nrow(toUpdate), " file(s)"))
print(paste0("Removed ", nrow(toRemove), " file(s)"))
