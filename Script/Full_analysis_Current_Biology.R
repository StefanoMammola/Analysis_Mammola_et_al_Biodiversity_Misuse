## ------------------------------------------------------------------------
## 'How much biodiversity is concealed in the word “biodiversity”?'
## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
# 'R script to reproduce the full analysis'
## ------------------------------------------------------------------------

# Analysis performed with R (v. R 4.1.0) and R studio (v. 1.4.1103)
# Authors: Stefano Mammola

####################
# Data preparation #
####################

# Clean the workspace -----------------------------------------------------

rm(list = ls())

# Loading R package -------------------------------------------------------

library("Amelia")
library("dplyr")
library("gam")
library("ggplot2")
library("MuMIn")
library("parameters")
library("performance")
library("psych")
library("PupillometryR")
library("tidyr")
library("tidytext")
library("tidyverse")
library("wordcloud")
library("RColorBrewer")

# Source functions and plot parameters ------------------------------------

source("Script/Functions.R")

# Loading database --------------------------------------------------------

db <- read.csv(file = "Data/Database_Full_study.csv", sep = '\t', dec='.', header = TRUE, as.is = FALSE)

str(db)
dim(db)

# Selecting paper to analyse
db <- db %>% filter(Analysis == "yes") ; db <- droplevels(db)

###########################
# Data preparation part 1 #
###########################

# Checking factor levels
levels(db$Geography)
levels(db$Domain)
levels(db$Method_data_collection)

# Calculating proportion of biodiversity
db$Biodiversity_prop <- rowSums(db[,36:91]) / length(36:91)
db$Animals_prop      <- rowSums(db[,c(36:70,89)]) / length(c(36:70,89))
db$Plants_prop       <- rowSums(db[,c(71:75,90)]) / length(c(71:75,90))
db$Fungi_prop        <- rowSums(db[,c(76:84,91)]) / length(c(76:84,91))
db$Micro_prop        <- rowSums(db[,c(85:88)])    / length(c(85:88))

#Title descriptors
db$Title_fanciness <- rowSums(db[,c(23,24)])
table(db$Title_fanciness) #too few obs

# Calculating total number of specifics to the title
db$Title_adjecties <- as.factor(rowSums(db[,27:29]))
table(db$Title_adjecties)

# Number of biodiversity facets
db$Facets_biodiversity <- as.factor(rowSums(db[,31:34]))

# Facets
table(db[db$Facets_biodiversity != 0,]$Facets_biodiversity)/nrow(db[db$Facets_biodiversity != 0,])
table(db[db$Facets_biodiversity != 0,]$Taxonomic_div)/nrow(db[db$Facets_biodiversity != 0,])
table(db[db$Facets_biodiversity != 0,]$Phylogenetic_div)/nrow(db[db$Facets_biodiversity != 0,])
table(db[db$Facets_biodiversity != 0,]$Functional_div)/nrow(db[db$Facets_biodiversity != 0,])
table(db[db$Facets_biodiversity != 0,]$Other_div)/nrow(db[db$Facets_biodiversity != 0,])

# Missing data 
Amelia::missmap(db)

# Chcking N° of coutries in the reference list
db$author_country <- as.character(db$author_country)

country_diversity <- c()

for(i in 1:nrow(db)){
  
  diversity_i <- strsplit(db[i,]$author_country, split = ";")[[1]]
  diversity_i <- trimws(diversity_i, which = c("both"))
  diversity_i <- unique(diversity_i)
  
  country_diversity <- append(country_diversity, length(diversity_i))}

db <- cbind(db,country_diversity) 

###########################
# Data preparation part 2 #
###########################

# Correcting citations by year of publication ----------------------

# What is the trend of citation over time?
db %>% ggplot(aes(x = Publication_year, y = tot_cites)) + 
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "year of publication", y = "citation")+
  theme_custom()

# Modeling temporal trend
M0 <- gam::gam(tot_cites ~ s(Publication_year), family = poisson, data = db) #constraining dof to 2 so that you don't get negative predicted citations
summary(M0)

par(mar= c(rep(2,4)))
plot(M0, se = TRUE)

# Taking the residual of citations
db <- data.frame(db, citation_residuals = resid(M0, type="pearson"))

# Checking how good it predict by comparing with a simple normalization by year of publication
db %>%  ggplot(aes(x = tot_cites/Publication_year, y = citation_residuals)) + #normalize by year
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "citation / year", y = "citation residuals")+
  theme_custom()

# Correcting altmetrics by year of publication ----------------------------

db_alt <- db %>% select(year = Publication_year,
                        Altmetrics = Altmetrics,
                        row_ID) %>% na.omit
                 
db_alt %>% ggplot(aes(x = year, y = Altmetrics)) + 
  geom_point(size = 1, alpha = 0.7, color = "grey40")+
  labs(x = "year of publication", y = "altmetric score")+
  theme_custom()

# Modeling temporal trend
M1 <- gam::gam(Altmetrics ~ s(year), family = poisson, data = db_alt) #constraining dof to 2 so that you don't get negative predicted citations
summary(M1)

dev.off()
par(mar= c(rep(2,4)))
plot(M1, se = TRUE)

# Taking the residual of citations & merge with the original db
db_alt <- data.frame(db_alt, Altmetrics_residuals = resid(M1,type="pearson")) # Taking the residual of citations

db_alt <- db_alt %>% select(row_ID,Altmetrics_residuals)

db <- db %>% dplyr::left_join(db_alt, by = "row_ID") ; rm(db_alt, M0, M1)

###########################
# Data preparation part 3 #
###########################

# Assigning Impact factor for each publication ----------------------------

# Uploading the database with all Impact factors between 1997 and 2020 (from Clarivate analytics Journal Citation report)
my_files   <- paste("Data/IF/", list.files("Data/IF/"), sep = '')
all_IF     <- lapply(my_files, read.csv, header = TRUE, sep = ",") 

# watch out that 2019 and 2020 are new the new WSC format!!
str(all_IF[[23]])  #Wrong!

# Reading thhem manually:
all_IF[[23]] <- read.csv(file = my_files[23], header = TRUE, sep = "\t")
all_IF[[24]] <- read.csv(file = my_files[24], header = TRUE, sep = "\t")

str(all_IF[[23]]) #ok!
str(all_IF[[24]]) #ok!

# Generating a database list all journal and their impact factor for each year
my_IF <- list() ; year_IF <- 1997:2020

for (k in 1 : length(all_IF)){
  
  if (k < 23) {

  df <- all_IF[[k]]
  df <- df[,c(2,4)]
  colnames(df) <- c("journal", "IF")
  df[df == "Not Available"] <- "NA"
  
  df$IF         <- as.numeric(as.character(df$IF))

  df$journal    <- sapply(df$journal, word.cleaner, remove.punctuation = TRUE) #simplify journal string
  df$journal    <- paste(df$journal, rep(year_IF[[k]],nrow(df)),sep='') #paste journal and year
  
  df <- df %>% distinct()
  
  my_IF[[k]] <- df
  
  } else {  
    df <- all_IF[[k]]

    colnames(df) <- c("journal", "IF")
    df[df == "Not Available"] <- "NA"
    
    df$IF         <- as.numeric(as.character(df$IF))
    df$journal    <- as.character(df$journal)
    
    df$journal    <- sapply(df$journal, word.cleaner, remove.punctuation = TRUE) #simplify journal string
    df$journal    <- paste(df$journal, rep(year_IF[[k]],nrow(df)),sep='') #paste journal and year
  
    df <- df %>% distinct()
    
    my_IF[[k]] <- df
  }
  
} 
#warnings() # ( Warning is due to creation of missing data)

# Unlist the list at last
IF_1997_2020 <- do.call("rbind", my_IF)
colnames(IF_1997_2020)[1] <- "JI_PY"

db$JI_PY <- paste(word.cleaner(db$journal,remove.punctuation = TRUE), db$Publication_year,sep='')

db <- db %>% dplyr::left_join(IF_1997_2020, by = "JI_PY")

rm(my_files,all_IF,my_IF,df,IF_1997_2020,k) #cleaning

##################
# Summary stats  #
##################

# How many studies have no biodiversity proportion?
nrow(db[db$Biodiversity_prop == 0,])/nrow(db) * 100 #22.1% do not consider any biodiversity group

# Removing study with no Biodiversity 
db2 <- db %>% filter(Biodiversity_prop > 0)

#####################################################
# Data preparation: part 2 (after data exploration) #
#####################################################

# Converting multiples levels in factors
method_split <- strsplit(as.character(db2$Method_data_collection), ";")

method <- c()
for(i in 1:length(method_split))
  method <- c(method, ifelse(length(method_split[[i]]) > 1, "Multiple", method_split[[i]]) )

geography_split <- strsplit(as.character(db2$Geography), ";")

geography <- c()
for(i in 1:length(geography_split))
  geography <- c(geography, ifelse(length(geography_split[[i]]) > 1, "Global", geography_split[[i]]) )

domain_split <- strsplit(as.character(db2$Domain), ";")

domain <- c()
for(i in 1:length(domain_split))
  domain <- c(domain, ifelse(length(domain_split[[i]]) > 1, "Multiple", domain_split[[i]]) )

db2$Method    <- method
db2$Geography <- geography
db2$Domain    <- domain

db2 <- db2 %>% mutate_at(vars("Method","Geography","Domain","Title_geo","Title_hab","Title_taxon"), as_factor)

# Set baseline
db2 <- within(db2, Geography <- relevel(Geography, ref = "Global"))
db2 <- within(db2, Domain    <- relevel(Domain,    ref = "Multiple"))
db2 <- within(db2, Method    <- relevel(Method,    ref = "Multiple"))

#clean
rm(method_split, method, geography, geography_split, domain, domain_split, i)

####################
# Summary stats #2 #
####################

# Biodiversity facets (%)
round((table(db2$Facets_biodiversity)/sum(table(db2$Facets_biodiversity)))*100,2)
table(db2$Phylogenetic_div)
table(db2$Functional_div)
table(db2$Taxonomic_div)
table(db2$Other_div)

#What proportion of biodiversity across studies?
mean(db2$Biodiversity_prop, na.rm = TRUE) * 100 #mean 
std(db2$Biodiversity_prop)*100 #sW 
range(db2$Biodiversity_prop, na.rm = TRUE) * 100 #range
getmode(db2$Biodiversity_prop)*100 #mode

# Checking temporal distribution

Percentile <- ifelse(db2$Biodiversity_prop > quantile(db2$Biodiversity_prop, c(.75)),"75–100 percentile","0–75 percentile")

model_0 <- glm(Biodiversity_prop ~ Publication_year, data = db2, family = quasibinomial(link = "logit"))
model_1 <- glm(Biodiversity_prop ~ Publication_year, data = db2[Percentile == "75–100 percentile", ], family = quasibinomial(link = "logit"))

summary(model_0)
summary(model_1)

Col_custom <- c(rev(RColorBrewer::brewer.pal(5, "Blues"))[1], "orange")

(plot1a <- ggplot(data = db2, aes(x = Publication_year, y = Biodiversity_prop)) + 
    geom_point(size = 3, shape = 21, alpha = 0.25,
               col = "grey30",
               fill = "grey30") +
    
    
    geom_smooth(method = "glm", formula = y ~ x, col = rev(RColorBrewer::brewer.pal(5, "Blues"))[1],
                fill = rev(RColorBrewer::brewer.pal(5, "Blues"))[1], alpha = 0.5,
                method.args = list(family = quasibinomial(link = "logit"))) +
    
    geom_smooth(data = db2[Percentile == "75–100 percentile", ], aes(x = Publication_year, y = Biodiversity_prop),
                method = "glm", formula = y ~ x, col = rev(RColorBrewer::brewer.pal(5, "Blues"))[1],
                fill = rev(RColorBrewer::brewer.pal(5, "Blues"))[1], alpha = 0.25,linetype = "dashed",
                method.args = list(family = quasibinomial(link = "logit"))) +
    labs(x = NULL, y = Y.label)+
    theme_custom()+ theme(legend.position = c(0.25, 0.85),
                          legend.text.align = 0,
                          legend.title = element_text(size = 12)))

########################################################################
# Regression model ------------------------------------------
########################################################################

# Data exploration --------------------------------------------------------

# Subset
db_glm <- db2 %>% select(year = Publication_year,
                        n_aut,
                        citation_residuals,
                        Method,
                        Phylogenetic_div,      
                        Functional_div,
                        Other_div,
                        Title_geo,
                        Title_hab,
                        Title_taxon,
                        Title_adjecties,
                        Biodiversity_prop,
                        Geography,
                        Domain,
                        country_diversity)

# Checking balancing of factors
table(db_glm$Method) #Citizen science/simulation too few records
levels(db_glm$Method)[c(6,7)] <- "Other" 

table(db_glm$Geography) #Antartica a bit unbalanced

table(db_glm$Domain) #OK
levels(db_glm$Domain)[c(3,4)] <- "Aquatic" 

db_glm <- droplevels(db_glm)

# Dependent var
db_glm$prop  <- rowSums(db2[,36:91])
db_glm$total <- length(36:91)

db_outliers <- db2[db2$Biodiversity_prop > 0.3,]

quantile(db_outliers$Biodiversity_prop)
#22 and 25 phyla/division

db_outliers$title
db_outliers$Biodiversity_prop

#db_glm <- db_glm[db_glm$citation_residuals < 200,]
db_glm <- db_glm[db_glm$prop < 20,]

# Scale continuous variables
db_glm <- db_glm %>% mutate_at(vars(year,citation_residuals,country_diversity),scale)

# Modelling  --------------------------------------------------------------

# Set formula
model_1 <- as.formula("prop ~ 
                      year + 
                      Domain + Geography + Method +
                      Phylogenetic_div + Functional_div + Other_div +
                      Title_geo + Title_hab + Title_taxon")

#Initial model
m1  <- glm(model_1, data = db_glm, family = "poisson")
performance::check_overdispersion(m1) #overdispersed

#Refit with quasibinomial due to overdispersion 
m1b  <- MASS::glm.nb(model_1, data = db_glm)

parameters::parameters(m1b, df_method = "wald")
performance::check_collinearity(m1b)
performance::r2(m1b)

# PostHoc
pairs(emmeans::emmeans(m1b, ~ Domain), simple = c("Domain"))
pairs(emmeans::emmeans(m1b, ~ Geography), simple = "Geography")
pairs(emmeans::emmeans(m1b, ~ Method), simple = "Method")

# Estract estimates
Estimates_m1 <- 
  m1b %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

# Set variable order and rename
order_var1 <- c("Year of publication",
               "Domain [Terrestrial]",
               "Domain [Aquatic]",
               "Biogeography [Palearctic]",
               "Biogeography [Afrotropical]",
               "Biogeography [Indomalayan]",
               "Biogeography [Neartic]",
               "Biogeography [Australasian]",
               "Biogeography [Neotropical]",
               "Biogeography [Antartic]",
               "Method [Review/Opinion]",
               "Method [Field sampling]",
               "Method [Big data]",
               "Method [Other]",
               "Phylogenetic diversity [yes]",
               "Functional diversity [yes]",
               "Other diversity [yes]",
               "Mention of location in title [yes]",
               "Mention of habitat in title [yes]",
               "Mention of taxon/a in title [yes]")

Estimates_m1$Variable <- order_var1 #Rename
Estimates_m1$Variable <- factor(Estimates_m1$Variable, rev(order_var1)) #Sort

sign <- ifelse(Estimates_m1$p > 0.05, "", ifelse(Estimates_m1$Estimate>0.01," *", " **")) #Significance
col_p <- ifelse(Estimates_m1$p > 0.05, "grey5", ifelse(Estimates_m1$Estimate>0,"orange",Col_custom[1])) #Significance

# Plot
plot1b <- ggplot2::ggplot(data = Estimates_m1) +

    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate),
              label = paste0(round(Estimates_m1$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
   
    labs(title = paste0("Dependent variable: Sampled biodiversity [N = ",nrow(db_glm),"]"),
         y = expression(paste("Effect size" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()+
  annotate(geom = 'text', x = 2, y = -0.9, size =5,
           label = paste0("R^2 ==",round(as.numeric(performance::r2(m1b)[1]),2)), parse = TRUE)

#############################
# Modelling article impact ##
#############################

#missing impact factor

# Preparing the data ------------------------------------------------------

db2$bio <- db2$Biodiversity_prop*56

db_impact <- db2 %>% select(year = Publication_year,
                        Altmetrics,
                        cit = tot_cites,
                        Altmetrics_residuals,
                        citation_residuals,
                        Biodiversity = bio,
                        Descriptors = Title_adjecties,
                        Title_geo,
                        Title_taxon,              
                        Title_hab,
                        n_aut,
                        journal,
                        IF,
                        country_diversity) %>% 
                        mutate_at(vars(starts_with("Title_"), journal), as_factor) 
#%>% 
#                        mutate_at(vars(n_aut, year, country_diversity, IF), scale)

# Setting baseline
levels(db_impact$Title_taxon) <- c("Mention","No mention")
levels(db_impact$Title_hab)   <- c("Mention","No mention")
levels(db_impact$Title_geo)   <- c("Mention","No mention")
levels(db_impact$Descriptors)  <- c("None","One descriptor", "Two descriptors", "Three descriptors")

db_impact <- within(db_impact, Title_geo   <- relevel(Title_geo, ref = "Mention"))
db_impact <- within(db_impact, Title_taxon <- relevel(Title_taxon, ref = "Mention"))
db_impact <- within(db_impact, Title_hab   <- relevel(Title_hab, ref = "Mention"))

#Random factor levels
nlevels(db_impact$journal) #281

# Collinearity
psych::pairs.panels(db_impact %>% as.data.frame %>% select(n_aut, country_diversity, IF, Altmetrics_residuals, citation_residuals))
#Dropping country diversity

# Checking outliers
#dotchart(db_impact$Biodiversity) #1 outliers
#dotchart(db_impact$country_diversity)
#dotchart(db_impact$n_aut)
#dotchart(db_impact$IF) 
db_impact$Biodiversity      <- log(db_impact$Biodiversity+1) #log transform
db_impact$IFlog             <- log(db_impact$IF+1) #log transform
db_impact$n_aut             <- log(db_impact$n_aut+1) #log transform
db_impact$country_diversity <- log(db_impact$country_diversity+1) #log transform

# What's the relationship between altmetrics and citations
db_impact %>% ggplot2::ggplot(aes(x = log(Altmetrics_residuals+1), y = log(citation_residuals+1))) + 
  geom_point(col = "grey10", fill = "grey30", size = 5, shape = 21, alpha = 0.3) + 
  geom_smooth(method = "lm",  se = TRUE, col = "blue", fill = "blue",
              formula = y ~ x) +
  theme_custom()

# Testing Citations -----------------------------------------------------------

#Check outliers
par(mar= c(rep(2,4)))
#dotchart(db_impact$citation_residuals)  #1 outlier

db_cit <- db_impact %>% select(citation_residuals,
                            Biodiversity,
                            Descriptors,
                            Title_geo,
                            Title_taxon,              
                            Title_hab,
                            country_diversity,
                            journal,
                            IFlog) %>% filter(citation_residuals < 150)

db_cit <- na.omit(db_cit)

# First model (all descriptors sum,ed)
model_3 <- as.formula("citation_residuals ~ Title_geo + 
                                            Title_taxon + 
                                            Title_hab + 
                                            country_diversity +
                                            IFlog +
                                            Biodiversity : Descriptors") #(1|journal)"

m3  <- lm(model_3, data = db_impact)
parameters::parameters(m3)
performance::r2(m3)

performance::check_model(m3)

#normality or residuals
hist(resid(m3), breaks = 200) #ok, just a bit of tail

#Contrast
pairs(emmeans::emmeans(m3, ~ Biodiversity : Descriptors), simple=c("Descriptors"))

# Estract estimates
Estimates_m3 <- 
  m3 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, t = 4) #rename

# Set variable order and rename
order_var3 <- c("Mention of location in title [no]",
                "Mention of taxon/a in title [no]",
                "Mention of habitat in title [no]",
                "Number of coauthors' countries",
                "Impact factor",
                "Sampled biodiversity : No descriptors",
                "Sampled biodiversity : 1 descriptor",
                "Sampled biodiversity : 2 descriptors", 
                "Sampled biodiversity : 3 descriptors")

Estimates_m3$Variable <- order_var3 #Rename
Estimates_m3$Variable <- factor(Estimates_m3$Variable, rev(order_var3)) #Sort

par <- parameters::parameters(m3) %>%
  dplyr::filter(!row_number() %in% 1)

par <- na.omit(par$p)

sign <- ifelse(par > 0.05, "", ifelse(par > 0.01," *", " **")) #Significance
col_p <- ifelse(par > 0.05, "grey5", Col_custom[1])

#Set variable order and rename

(plot1c <- ggplot2::ggplot(data = Estimates_m3) +
    
    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_text(aes(Variable, Estimate), label = paste0(round(Estimates_m3$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
    
    labs(title = paste0("Dependent variable: Citations [N = ", nrow(db_cit),"]"),
         y = expression(paste("Effect size" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip() +
    annotate(geom = 'text', x = 1, y = 7.5, size = 5,
             label = paste0("R^2 ==",round(as.numeric(performance::r2(m3)[1]),2)), parse = TRUE)
)

# Check the interaction
(plot_interaction1 <- db_cit %>% ggplot2::ggplot(aes(x = Biodiversity, y = citation_residuals)) + 
    facet_wrap( ~ Descriptors, nrow = 2, ncol = 2) +
    geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
    geom_smooth(method = "lm",  se = TRUE, col = Col_custom[1], fill = Col_custom[1],
                formula = y ~ x) +
    labs(x = "Sampled biodiversity [log-transformed]", 
         y = "Citations [residuals]",
         title = "Interaction sampled biodiversity * N° of descriptors")+ theme_custom())

# Testing Altmetric -------------------------------------------------------

db_alt <- db_impact %>% select(Altmetrics_residuals,
                            Biodiversity,
                            Descriptors,
                            Title_geo,
                            Title_taxon,              
                            Title_hab,
                            country_diversity,
                            journal,
                            IFlog)

db_alt <- na.omit(db_alt)

db_alt$Descriptors <- droplevels(db_alt$Descriptors)

# Checking outliers
# par(mar= c(rep(2,4)))
# dotchart(db_alt$Altmetrics_residuals) # OK

# Set formula
model_4 <- as.formula("Altmetrics_residuals ~ Title_geo + 
                                              Title_taxon + 
                                              Title_hab + 
                                              country_diversity +
                                              IFlog +
                                              Biodiversity : Descriptors ")


# Initial model
m4 <- lm(model_4, data = db_alt)
performance::r2(m4)
parameters::parameters(m4)

#Collinearity
performance::check_model(m4)

#normality or residuals
hist(resid(m4), breaks = 200) #ok, just a bit of tail

#Contrast
pairs(emmeans::emmeans(m4, ~ Biodiversity : Descriptors), simple=c("Descriptors"))

# Estract estimates
Estimates_m4 <- 
  m4 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, t = 4) #rename

# Set variable order and rename
order_var4 <- c("Mention of location in title [no]",
                "Mention of taxon/a in title [no]",
                "Mention of habitat in title [no]",
                "Number of coauthors' countries",
                "Impact factor",
                "Sampled biodiversity : No descriptors",
                "Sampled biodiversity : 1 descriptor",
                "Sampled biodiversity : 2 descriptors", 
                "Sampled biodiversity : 3 descriptors")

Estimates_m4$Variable <- order_var4 #Rename
Estimates_m4$Variable <- factor(Estimates_m4$Variable, rev(order_var4)) #Sort

par <- parameters::parameters(m4) %>%
  dplyr::filter(!row_number() %in% 1)

par <- na.omit(par$p)

sign <- ifelse(par > 0.05, "", ifelse(par > 0.01," *", " **")) #Significance
col_p <- ifelse(par > 0.05, "grey5", Col_custom[1])

(plot1d <- ggplot2::ggplot(data = Estimates_m4) +
    
    geom_pointrange(aes(x = Variable, 
                        y = Estimate,
                        ymin = Estimate-SE, 
                        ymax = Estimate+SE), col = col_p, size = 0.5) + 
    
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    
    geom_text(aes(Variable, Estimate), label = paste0(round(Estimates_m4$Estimate,2),sign), 
              vjust = -1, size = 3, col = col_p) +
    
    labs(title = paste0("Dependent variable: Altmetric score [N = ", nrow(db_alt),"]"),
         y = expression(paste("Effect size" %+-% "Standard Error")),
         x = NULL)+
    theme_custom() + theme(axis.text.y  = element_text(colour = rev(col_p))) + coord_flip()+
    annotate(geom = 'text', x = 1, y = 7.5, size = 5,
             label = paste0("R^2 ==",round(as.numeric(performance::r2(m4)[1]),2)), parse = TRUE)
)

# Check the interaction
(plot_interaction2 <- db_alt %>% ggplot2::ggplot(aes(x = Biodiversity, y = Altmetrics_residuals)) + 
    facet_wrap( ~ Descriptors, nrow = 2, ncol = 2) +
    geom_point(col = "grey10", fill = "grey30", size = 3, shape = 21, alpha = 0.3)+
    geom_smooth(method = "lm",  se = TRUE, col = Col_custom[1], fill = Col_custom[1],
                formula = y ~ x) +
    labs(x = "Sampled biodiversity [log-transformed]", 
         y = "Altmetric score [residuals]",
         title = NULL)+ theme_custom())

# Store the final figure
pdf(file = "Figure_Current_Biology/Figure_1.pdf", width = 18, height = 15)

ggpubr::ggarrange(plot1a,plot1b,plot1c,plot1d,
                  
                  common.legend = FALSE,
                  align = "h",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2) #warnings() are due to NA removal
dev.off()




# Supplementary figure S2 ----------------------------------------------------------------

vector <- sort(apply(db2[,36:86],2, sum, na.rm = TRUE), decreasing = TRUE)

bar1 <- data.frame(Phyla = names(vector), N = vector) 

threshold <- 10 #cut of to merge in multiple category
Taxa_to_rename <- ifelse(bar1$N < threshold, bar1$N, NA)
N_to_rename    <- sum(Taxa_to_rename,na.rm = TRUE)
Taxa_to_rename <- length(na.omit(Taxa_to_rename))

bar1 <- bar1[ 1 : (nrow(bar1)-Taxa_to_rename), ]

bar1 <- rbind(bar1, data.frame(Phyla = paste0("Others (n = ",Taxa_to_rename,")"), N = N_to_rename))

bar1 <- cbind(bar1, Type = c(rep("Animals",2),
                             "Plants",
                             rep("Animals",2),
                             "Plants",
                             "Microorganisms",
                             "Plants",
                             rep("Animals",3),
                             rep("Plants",2),
                             "Animals",
                             "Fungi", #ascoo
                             rep("Animals",2),
                             "Multiple",
                             "Fungi",
                             "Multiple"))

bar1$Phyla <- factor(bar1$Phyla,levels = bar1$Phyla)
levels(bar1$Phyla)[18] <- "Protista *"

col_bar <- rev(RColorBrewer::brewer.pal(5, "Blues"))

barS1 <- ggplot(bar1, aes(x = Phyla, y = N))+
  geom_bar(aes(fill= Type),stat="identity", alpha=1, colour = "black")  +
  labs(x = NULL, y = "Count")+
  scale_fill_manual(values = col_bar)+
  theme_custom()+ 
  theme(legend.position = c(0.75,0.45))+ 
  coord_flip()+
  annotate(geom = "text", y = bar1[nrow(bar1),2] + 2, 
           x = nrow(bar1), 
           label = paste0("Phyla/division with sample size < ", threshold),
           hjust = 0)

rm(Taxa_to_rename, N_to_rename, vector, bar1, threshold)

## Only with titles containong no descriptors

bar2 <- db2[db2$Title_adjecties == 0,]

vector <- sort(apply(bar2[,36:86],2, sum, na.rm = TRUE), decreasing = TRUE)

bar2 <- data.frame(Phyla = names(vector), N = vector)

threshold2 <- 3 #cut of to merge in multiple category
Taxa_to_rename <- ifelse(bar2$N < threshold2, bar2$N, NA)
N_to_rename    <- sum(Taxa_to_rename,na.rm=TRUE)
Taxa_to_rename <- length(na.omit(Taxa_to_rename))

bar2 <- bar2[ 1 : (nrow(bar2)-Taxa_to_rename), ]

bar2 <- rbind(bar2, data.frame(Phyla = paste0("Others (n = ",Taxa_to_rename,")"), N = N_to_rename))

bar2 <- cbind(bar2, Type = c(rep("Animals",2),
                             rep("Plants",2),
                             "Microorganisms",
                             rep("Animals",3),
                             rep("Plants",2),
                             rep("Animals",2),
                             "Multiple"))

bar2$Phyla <- factor(bar2$Phyla, levels = bar2$Phyla)

col_bar2 <- col_bar[c(1,3:5)]

barS2 <- ggplot(bar2, aes(x = Phyla, y = N))+
  geom_bar(aes(fill = Type),stat="identity", alpha=1, colour = "black")  +
  labs(x = NULL, y = "Count")+
  scale_fill_manual(values = col_bar2) +
  theme_custom() + 
  theme(legend.position = "none") + 
  coord_flip() +
  annotate(geom = "text", y = bar2[nrow(bar2),2] + 1, 
           x = nrow(bar2), 
           label = paste0("Phyla/division with sample size < ",threshold2),
           hjust = 0)

# Checking distribution by region, method and system
(plotS1a <- db2 %>% 
    drop_na(Geography,Biodiversity_prop) %>% 
    #filter(Biodiversity_prop < 0.4) %>%  #removing 1 outlier
    ggplot(aes(x = Geography, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                     alpha = 1, fill = RColorBrewer::brewer.pal(5, "Blues")[5], col = "white", adjust = 1.5) +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = RColorBrewer::brewer.pal(5, "Blues")[3]) +
    geom_boxplot(width = 0.2,  col = RColorBrewer::brewer.pal(5, "Blues")[5], outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom() + theme(axis.text.y = element_text(size = 15)) + coord_flip())

(plotS1b <- db2 %>% 
    drop_na(Domain,Biodiversity_prop) %>% 
    ggplot(aes(x = Domain, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                     alpha = 1, fill= RColorBrewer::brewer.pal(5, "Blues")[5], col = "white", adjust = 1.5) +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = RColorBrewer::brewer.pal(5, "Blues")[3]) +
    geom_boxplot(width = 0.2, col = RColorBrewer::brewer.pal(5, "Blues")[5], outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom() + theme(axis.text.y = element_text(size = 15)) + coord_flip())

levels(db2$Method)[c(6,7)] <- "Other"

(plotS1c <- db2 %>% 
    drop_na(Method,Biodiversity_prop) %>% 
    ggplot(aes(x = Method, y = Biodiversity_prop)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                     alpha = 1, fill= RColorBrewer::brewer.pal(5, "Blues")[5], col = "white", adjust = 1.5) +
    geom_point(position = position_jitter(width = 0.15), size = 1, alpha = 0.7, color = RColorBrewer::brewer.pal(5, "Blues")[3]) +
    geom_boxplot(width = 0.2, col = RColorBrewer::brewer.pal(5, "Blues")[5], outlier.shape = NA, alpha = 0) +
    labs(y = Y.label, x = NULL) +
    theme_custom() + theme(axis.text.y = element_text(size = 15)) + coord_flip())

pdf(file = "Figure_Current_Biology/Figure_S2.pdf", width = 14, height = 10)

ggpubr::ggarrange(plotS1a, plotS1b, plotS1c, barS1,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2) #warnings() are due to NA removal

dev.off()
#end