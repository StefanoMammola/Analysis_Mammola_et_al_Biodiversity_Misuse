## ------------------------------------------------------------------------
## 'Script to standardize species distribution (by continent)
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.0.3) and R studio (v. 1.4.1103)

# Loading R packages ------------------------------------------------------

library("arakno")
library("dplyr")
library("countrycode")
library("tidyverse")

# Loading the database ----------------------------------------------------

db <- read.csv(file = "Data/db_initial.csv", header = TRUE, sep = "\t", as.is = FALSE)

# Extract ISO codes -------------------------------------------------------

# ISO code for the World Spider Catalog
ISO <- force(wscmap)
colnames(ISO)[1] <- "distribution"

# ISO code vs country
data_cont <- force(codelist)
data_cont <- data_cont %>% select(continent,code=genc3c)
data_cont <- na.omit(data_cont)

Cont <- unique(data_cont$continent) #Continent string

# Match ISO code vs distribution in the main database
db$distribution <- as.character(db$distribution)

distribution_string <- lapply(1:nrow(db), function (x) tolower(strsplit(db$distribution[x], 
                       c("\\, | and | to |from |\\?|\\/|probably|possibly|introduced"))[[1]]) )

# Takes ca. 2 minute:
list_distr <- lapply(1:length(distribution_string), 
                     function (x) colSums(ISO[ISO$distribution %in% distribution_string[[x]] , ][,2:ncol(ISO)]))

distribution <- do.call("rbind",list_distr)
head(distribution)

# Associate continent
continent_list <- list()

for(i in 1:nrow(distribution))  {
  distr_i <- distribution[i, ]
  continent <- unique(data_cont[data_cont$code %in% names(distr_i[distr_i > 0]),]$continent) 
  continent <- factor(continent, levels = Cont)
  continent_list[[i]] <- table(continent)
  }

# Store result
final_distribution <- data.frame(distribution = db$distribution, 
                                 do.call("rbind",continent_list))

head(final_distribution)

#Attach to the main database
db <- dplyr::left_join(db,final_distribution, by = "distirbution")

# (Note: Afterward, the final assignment have been manually checked to correct mistakes and unmatched distributions)