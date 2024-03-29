## ------------------------------------------------------------------------
## ''How much biodiversity is concealed in the word “biodiversity”?'
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola

# Loading R package -------------------------------------------------------

library("xlsx")
library("wosr")

# WoS queries -------------------------------------------------------------

# Data source: Web of Science, accessed on 22.04.2021
# [Helsinki, Google Chrome, macOS High Sierra 10.13.6]

#Setting sid
sid <- auth(NULL, password = NULL) #change with your WoS access

#Setting WoS collection
coll <-
  c("SCI", "SSCI", "AHCI", "ISTP", "ISSHP", "BSCI", "BHCI", "ESCI")

# WoS query ---------------------------------------------------------------

query <-
  paste(
    '(TI = "Biodiversity") AND (DT = "Article") AND PY = (1965-2020) AND (WC = "Ecology" OR WC = "Soil Science" OR WC = "Environmental Studies" OR WC = "Environmental Sciences"
                         OR WC = "Marine & Freshwater Biology" OR WC = "Multidisciplinary Sciences" OR WC = "Paleontology")'
  )

(n_tot <- query_wos((query), editions = coll, sid = sid))
download_1 <- pull_wos(query, editions = coll, sid = sid)

# Organizing the data -----------------------------------------------------

pub      <- download_1$publication
aut      <- download_1$author
aff      <- download_1$address
key      <- download_1$keyword
key1     <- download_1$keywords_plus
WC_cat   <- download_1$jsc
doc_type <- download_1$doc_type

authors        <- c()
n_aut          <- c() 
authors_id     <- c()
author_country <- c()
keywords       <- c()
keywords_plus  <- c()
WC_category    <- c()
type           <- c()

for (i in 1:nrow(pub)){
  
  #extracting the id
  id <- pub[i,]$ut
  
  #storing authors
  aut_i <- aut[aut$ut == id,]
  authors <- c(authors, paste(aut_i$display_name, collapse =" ; "))
  
  n_aut   <- append(n_aut, nrow(aut_i))
  authors_id <- c(authors_id, paste(aut_i$daisng_id, collapse = " ; "))
  
  #store affiliation
  aff_i <- aff[aff$ut == id,]
  author_country <- c(author_country, paste(aff_i$country, collapse = " ; "))
  
  #keywords
  key_i <- key[key$ut == id,]
  keywords <- c(keywords, paste(key_i$keyword, collapse = " ; "))
  
}

final_db <- data.frame(pub,
                       authors,
                       n_aut, 
                       authors_id,
                       author_country,
                       keywords,
                       year = format(as.POSIXct(pub$date, format = "%m/%d/%Y"), format = "%Y"))

# Saving it ---------------------------------------------------------------

xlsx::write.xlsx(final_db, "Data/Biodiversity_WOS_V0.xlsx")

# Random sample - Pilot study ---------------------------------------------
n_part  <- 15
range_year <- seq(from=2000,to=2020,by=1) #21 study each, between 2000 and 2020

names <- c("Fabio Cianferoni",
           "Diego Fontaneto",
           "Emiliano Mori",
           "Ilaria Rosati",
           "Simonen Tenan",
           "Carmelo Fruciano",
           "Paola Pollegioni",
           "Michelangelo Morganti",
           "Fernando Urbano",
           "Núria Macias",
           "Jagoba Malumbres-Olarte",
           "Ana Lozano",
           "Marija Miličić",
           "Dinarte Teixeira", #note he later left the study
           "Carol Fukushima")

index <- c()

for (i in 1:length(range_year)){
  db_year <- final_db[final_db$year == range_year[i], ]
  index   <- append(index, sample(db_year$ut)[1 : n_part])
}

db_sampled   <- final_db[final_db$ut %in% index,]
db_unsampled <- final_db[!final_db$ut %in% index,]

nrow(db_sampled) + nrow(db_unsampled) #check

#Creating a new database
db_sampled   <- data.frame(db_sampled,   Assignment = rep(names, length(range_year)))
db_unsampled <- data.frame(db_unsampled, Assignment = rep("unsampled",nrow(db_unsampled)))

final_db2 <- rbind(db_sampled,db_unsampled)

# Saving it ---------------------------------------------------------------

xlsx::write.xlsx(final_db2, "Data/Biodiversity_WOS_V1.xlsx")

# Random sampling - Full study --------------------------------------------

db_assigned  <- final_db2[final_db2$Assignment != "unsampled",]
db_to_assign <- final_db2[final_db2$Assignment == "unsampled",]

n_part  <- 19

names_part <- c("Fabio Cianferoni",
                "Diego Fontaneto",
                "Emiliano Mori",
                "Ilaria Rosati",
                "Simone Tenan",
                "Carmelo Fruciano",
                "Paola Pollegioni",
                "Michelangelo Morganti",
                "Fernando Urbano",
                "Núria Macias",
                "Jagoba Malumbres-Olarte",
                "Ana Munévar",
                "Marija Miličić",
                "Dinarte Teixeira", #note he later left the study
                "Carol Fukushima", 
                "Paolo Domenici",
                "Lucia Bongiorni",
                "Angelina Lo Giudice",
                "Gemma Biondo")

n_to_samp <- ceiling(sum((c(rep(30,15),rep(50,4))))/20) #50 study each

index <- c()

for (i in 1:length(range_year)){
  
  db_year <- db_to_assign[db_to_assign$year == range_year[i], ]
  index   <- append(index, sample(db_to_assign$ut)[1 : n_to_samp])
}

db_to_assign_sampled   <- db_to_assign[db_to_assign$ut %in% index,]
db_to_assign_unsampled <- db_to_assign[!db_to_assign$ut %in% index,]

nrow(db_to_assign_sampled) + nrow(db_to_assign_unsampled) ; nrow(db_to_assign) #check

db_to_assign_sampled$Assignment <- c(sample(c(rep(names_part,30),
                                              rep(names_part[16:19],20))), #people with 50 
                                     rep("unsampled",16))

final_db3 <- rbind(db_assigned,db_to_assign_sampled,db_to_assign_unsampled)

#check
nrow(final_db3) ; nrow(final_db)
sort(table(final_db3$Assignment))

xlsx::write.xlsx(final_db2, "Data/Biodiversity_WOS_V1.xlsx")
