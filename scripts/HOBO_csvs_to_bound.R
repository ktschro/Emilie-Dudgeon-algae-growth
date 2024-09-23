#script to read in HOBO files and add them all to one dataframe. This works for all HOBO files that have been saved as a csv with their file names as something descriptive to the experimental treatment

#set wd and read in packages
setwd("/Users/katieschroeder/Documents/GitHub/Emilie-Dudgeon-algae-growth")
library(tidyverse)

csv_names <- list.files(path = "HOBO-csvs/", #set the path to your folder with csv files
                        pattern = "*.csv", #select all csv files in the folder
                        full.names = T) #output full file names (with path)

#id for joining
csv_names2 <- data.frame(file = csv_names, 
                         id = as.character(1:length(csv_names))) 

hobo_all <- csv_names %>% 
  lapply(read_csv) %>% #read all the files at once
  bind_rows(.id = "id") %>% #bind all tables into one object, and give id for each
  left_join(csv_names2) #join month column created earlier

#edit the file name to get just the logger ID
#the specific strings used in gsub will vary depending on your file names
hobo_all$logger_ID <- gsub(".*HOBO-csvs//(.+)_emilie.csv","\\1",hobo_all$file)

write.csv(hobo_all,"edited-files/hobo_all.csv")
