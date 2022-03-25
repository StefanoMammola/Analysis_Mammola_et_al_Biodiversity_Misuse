###############################################################

## An expert-curated global database on online newspaper articles on spiders and spider bites
## Mammola, S. et al. 2021

###############################################################

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Last update: 29 Jul 2021, Helsinki, Finland
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

###############################################################

# clean the workspace -----------------------------------------------------

rm(list=ls())

# Working directory -------------------------------------------------------

setwd("/Users/stefanomammola/Desktop/PAPERS IN CORSO/SPIDER NEWS GLOBAL/Mammola et al DataPaper") #change with your working directory

# Loading R package -------------------------------------------------------

library("Amelia")        # A Program for Missing Data  
library("countrycode")   # Convert Country Names and Country Codes   
library("dplyr")         # A Grammar of Data Manipulation         
library("ggalt")         # Extra Coordinate Systems, 'Geoms', Statistical Transformations, Scales and Fonts for 'ggplot2'         
library("ggplot2")       # Create Elegant Data Visualisations Using the Grammar of Graphics       
library("ggthemes")      # Extra Themes, Scales and Geoms for 'ggplot2' 
library("grid")          # For the silhouhette
library("gridExtra")     # Miscellaneous Functions for "Grid" Graphics     
library("maps")          # Draw Geographical Maps
library("psych")         # For Cohen's kappa
library("png")           # For danger symbol
library("PupillometryR") # A Unified Pipeline for Pupillometry Data 
library("tidytext")      # Text Mining using 'dplyr', 'ggplot2', and Other Tidy Tools
library("tidyverse")     # Easily Install and Load the 'Tidyverse'     
library("wordcloud")     # Word Clouds

# Loading useful functions ------------------------------------------------

# Standard error
std <- function(x) sd(x, na.rm = TRUE) / sqrt(length(x))

# Plot style. 

#Modified from: https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/

theme_niwot <- function(){
    theme_bw() +
      theme(#text = element_text(family = "Arial"),
            axis.text = element_text(size = 14), 
            axis.title = element_text(size = 18),
            axis.line.x = element_line(color="black"), 
            axis.line.y = element_line(color="black"),
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),                                          
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),  
            plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
            plot.title = element_text(size = 18, vjust = 1, hjust = 0),
            legend.text = element_text(size = 12),          
            legend.title = element_blank(),                              
            legend.position = c(0.95, 0.15), 
            legend.key = element_blank(),
            legend.background = element_rect(color = "black", 
                                             fill = "transparent", 
                                             size = 2, linetype = "blank"))
  }

###############################################################

## Data preparation:

###############################################################

# Loading the Database ----------------------------------------------------

db <- read.csv(file = "Data_spider_news_global.csv", sep = '\t', dec = '.', header = TRUE, as.is = FALSE)

str(db)
dim(db)

#Converting factors to chr
db$Title <- as.character(db$Title)
db$Notes <- as.character(db$Notes)

#Sort by country of the search
db <- db %>% arrange(Country_search) 

# Creating new variables --------------------------------------------------

#Type of Event
db <- db %>% mutate(TypeEvent = Bite + Death)

db$TypeEvent <- as.factor(db$TypeEvent) ; levels(db$TypeEvent) <- c("Encounter","Bite","Deadly bite")

#Total number of Errors
db <- db %>% mutate(TotalError = tidyr::replace_na(as.numeric(Taxonomic_error),0) + 
                                 tidyr::replace_na(as.numeric(Venom_error),0)     + 
                                 tidyr::replace_na(as.numeric(Anatomy_error),0)  + 
                                 tidyr::replace_na(as.numeric(Photo_error),0))

#Expert
db <- db %>% mutate(TotalExpert = tidyr::replace_na(as.numeric(Expert_arachnologist),0) + 
                                  tidyr::replace_na(as.numeric(Expert_doctor),0)        + 
                                  tidyr::replace_na(as.numeric(Expert_others),0))

db$TotalExpert <- ifelse(db$TotalExpert > 0, 1, 0) #converting to binary

#Date
db$Year_news <- as.Date(paste(db$d,db$m, db$yr, sep = "/"), format = '%d/%m/%Y')

#Database only with distinct events
db_unique_event <- distinct(db, ID_Event, .keep_all = TRUE) 

#Database only with distinct news
db_unique_news <- distinct(db, ID, .keep_all = TRUE) 

###############################################################

## Summary statistics:

###############################################################

# Missing data -------------------------------------------------------------

Amelia::missmap(db)

# Summary statistics ------------------------------------------------------

#Distinct news
nrow(db_unique_news)

#Distinct event
nrow(db_unique_event)
table(db_unique_event$TypeEvent)

#Number of unique species involved
nlevels(droplevels(db$Species))

as.character(sort(unique(droplevels(db$Species))))

#Number of country and news country
nlevels(unique(droplevels(db_unique_news$Country_search)))
sort(table(db_unique_news$Country_search))

#Number of lenguages
nlevels(unique(droplevels(db$Lenguage))) 

# % of news with error
sum(ifelse(db_unique_news$TotalError > 0, 1, 0), na.rm = TRUE) / nrow(db_unique_news) * 100
sum(db_unique_news$Taxonomic_error, na.rm = TRUE) / nrow(db_unique_news) * 100
sum(db_unique_news$Venom_error, na.rm = TRUE) / nrow(db_unique_news) * 100
sum(db_unique_news$Anatomy_error, na.rm = TRUE) / nrow(db_unique_news) * 100
sum(as.numeric(db_unique_news$Photo_error), na.rm = TRUE) / nrow(db_unique_news) * 100

# % of sensationalistic news
sum(na.omit(db_unique_news$Sensationalism)) / nrow(db) * 100

# % of news with Expert
sum(db_unique_news$TotalExpert, na.rm = TRUE) / nrow(db) * 100
sum(db_unique_news$Expert_arachnologist, na.rm = TRUE) / nrow(db) * 100
sum(db_unique_news$Expert_doctor, na.rm = TRUE) / nrow(db) * 100
sum(db_unique_news$Expert_others, na.rm = TRUE) / nrow(db) * 100

# Mean latitude of bite and death
median(db[db$TypeEvent == "Encounter",]$lat, na.rm = TRUE) ; std(db[db$TypeEvent == "Encounter",]$lat)
median(db[db$TypeEvent == "Bite",]$lat, na.rm = TRUE) ; std(db[db$TypeEvent == "Bite",]$lat)
median(db[db$TypeEvent == "Deadly bite",]$lat, na.rm = TRUE) ; std(db[db$TypeEvent == "Deadly bite",]$lat)


###############################################################

## Data visualization:

###############################################################

# Temporal trends: Figure 2 -----------------------------------------------

# Panel a

# Annual trend
db_yr <- data.frame(table(db_unique_news$yr,db_unique_news$TypeEvent)) ; colnames(db_yr) <- c("yr","Event","N")

(plot_year <- ggplot(db_yr, aes(x=yr, y=as.numeric(N), fill=Event))+
    geom_bar(stat="identity",alpha=0.8,colour = "black")+
    scale_x_discrete(breaks = c(2010:2020), labels = as.character(2010:2020))+ 
    scale_fill_manual(values =  c("turquoise3", "orangered", "grey10"))+
    labs(title="(a)", 
         x=NULL, 
         y = "Number of unique news",
         subtitle = "Annual distribution of news")+
    annotate("segment", 
             x    = 11, 
             xend = 11, 
             y    = 720, 
             yend = 620, 
             colour = "black", 
             size = 0.5, 
             arrow=arrow(ends = "last",angle = 45, length = unit(.2,"cm")))+
    annotate("text", 
             x    = 11, 
             y    = 780, 
             label = "Partial\ndata",
             colour = "black")+
    theme_niwot()+
    theme(legend.position = c(0.1, 0.75),
          plot.subtitle = element_text(size = 12, color = "#222222"))
)

# Monthly trends
# Plot re-adapted from:
# https://github.com/juliacat23/tidytuesday/blob/main/Week%204%20-%20Netflix/Netflix.R

# Preparing dataset of distribution of news by month for the north and south hemisphere
db_month_N <- data.frame(table(db[db$lat3>0,]$m)) ; colnames(db_month_N) <- c("Month", "N")
db_month_N$Month <- as.factor(db_month_N$Month) ; levels(db_month_N$Month) <- month.abb

db_month_S <- data.frame(table(db[db$lat3<0,]$m)) ; colnames(db_month_S) <- c("Month", "N")
db_month_S$Month <- as.factor(db_month_S$Month) ; levels(db_month_S$Month) <- month.abb

# Panel b
(plot_month_N <- ggplot(db_month_N, aes(x = as.factor(Month), y = N)) +
    geom_col(fill = c(rep("grey60",6),rep("grey10",3),rep("grey60",3))) + 
    scale_fill_identity() + 
    coord_polar() +
    labs(x="", y="",
         title = "(b)", 
         subtitle = "Northern hemisphere [Latitude > 0°]",
         caption = NULL)+
    theme_minimal() + 
    theme( 
      plot.title = element_text(size = 18), 
      plot.title.position = "plot", 
      plot.subtitle = element_text(size = 12), 
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14))
)

# Panel c
(plot_month_S <- ggplot(db_month_S, aes(x = as.factor(Month), y = N)) +
    geom_col(fill = c("grey10",rep("grey60",2),rep("grey10",2),rep("grey60",4),rep("grey10",2),"grey60")) + 
    scale_fill_identity() + 
    coord_polar() +
    labs(x="", y="",
         title = "(c)", 
         subtitle = "Southern hemisphere [Latitude < 0°]",
         caption = NULL) +
    theme_minimal() + 
    theme( 
      plot.title = element_text(size = 18), 
      plot.title.position = "plot", 
      plot.subtitle = element_text(size = 12), 
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14))
)

#Arrange in a plate
pdf("Figure_2.pdf", width = 15, height = 7)

gridExtra::grid.arrange(plot_year,plot_month_N,plot_month_S, 
                       layout_matrix = rbind(c(1,1,2),
                                             c(1,1,3))
                       )

dev.off()

# Maps of the events & spider species: Figure 3 -------------------------------

# Panel a
world <- map_data("world")

(map1 <- ggplot() +
    geom_map(map = world, data = world,
             aes(long, lat, map_id = region), 
             color = "gray50", fill = "gray95", size = 0.3) +
    geom_point(data = db_unique_event %>% drop_na(TypeEvent), 
               aes(x = lon, y = lat,colour=as.factor(TypeEvent)),
               alpha = 0.8, size = 1,
               shape = 16)+
  ylim(-50,75)+
  labs(title = "(a)", fill = "Type of event: ", y = "Latitude [°]") +
  scale_colour_manual(values = c("turquoise3", "orangered", "grey10")) +
  theme_bw() +
    theme(
      axis.text.x  = element_text(size = 14,color="white"), 
      axis.text.y  = element_text(size = 14),
      axis.title.y = element_text(size = 18),
      axis.title.x = element_blank(), 
      axis.line.x = element_line(color="white"), 
      axis.line.y = element_line(color="black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),  
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 18, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = c(0.1, 0.2), 
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "white", 
                                       size = 2, linetype = "blank"))
)

# Latitudinal distribution - panel b
(plot_lat <- 
    ggplot(data = db_unique_event %>% drop_na(TypeEvent,lon), 
           aes(x = TypeEvent, y = lat, fill = TypeEvent)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
    geom_point(aes(y = lat, color = TypeEvent), 
    position = position_jitter(width = 0.15), size = 1, alpha = 0.2) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.4) +
    ylim(-50,75)+
    labs(title = "(b)", x = NULL) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values =  c("turquoise3", "orangered", "grey10")) +
    scale_colour_manual(values = c("turquoise3", "orangered", "grey10")) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(size = 14), 
      axis.text.y  = element_text(size = 14),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(), 
      axis.line.x = element_line(color="black"), 
      axis.line.y = element_line(color="black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),  
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 18, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = c(0.1, 0.2), 
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "black", 
                                       size = 2, linetype = "blank"))
  
)

# Distribution of news by family and genuses

# Preparing dataset of Families 
Bar_plot_family <- data.frame(sort(table(db$Family))) ; colnames(Bar_plot_family) <- c("Family","N")

#Summarize family with less thab 50 occurrence ina new category (others)
Vector <- as.character(Bar_plot_family$Family) ; Fam_to_rename <- sum(ifelse(Bar_plot_family$N>50,0,1))

renamed_Vector <- c(rep( paste("Others (n= ",Fam_to_rename,")",sep=''), Fam_to_rename), 
                    Vector[c(Fam_to_rename+1) : nrow(Bar_plot_family)]) 

Bar_plot_family <- data.frame(Bar_plot_family,Family_2 = as.factor(renamed_Vector));
Bar_plot_family$Family_2 <- factor(Bar_plot_family$Family_2,levels = renamed_Vector[Fam_to_rename : nrow(Bar_plot_family)]) ; rm(Vector,Fam_to_rename,renamed_Vector)
Bar_plot_family <- Bar_plot_family %>% group_by(Family_2) %>% summarise(N = sum(N))

# Preparing dataset of genuses
Bar_plot_genus <- data.frame(sort(table(db$Genus))) ; colnames(Bar_plot_genus) <- c("Genus","N")

#Summarize family with less than 20 occurrence ina new category (others)
Vector <- as.character(Bar_plot_genus$Genus) ; Gen_to_rename <- sum(ifelse(Bar_plot_genus$N>20,0,1))

renamed_Vector <- c(rep( paste("Others (n= ",Gen_to_rename,")",sep=''), Gen_to_rename), 
                    Vector[c(Gen_to_rename+1) : nrow(Bar_plot_genus)]) 

Bar_plot_genus <- data.frame(Bar_plot_genus, Genus_2 = as.factor(renamed_Vector))
Bar_plot_genus$Genus_2 <- factor(Bar_plot_genus$Genus_2,levels = renamed_Vector[Gen_to_rename : nrow(Bar_plot_genus)]) ; rm(Vector,Gen_to_rename,renamed_Vector)
Bar_plot_genus <- Bar_plot_genus %>% group_by(Genus_2) %>% summarise(N = sum(N))

# Panel c
col_fam = c(rep("grey50",nrow(Bar_plot_family)-4),"darkgoldenrod2","darkorchid4", "black", "brown3")

(plot_family <- ggplot(Bar_plot_family, aes(x= Family_2, y=N))+
    geom_bar(stat="identity",alpha=0.8,colour = "black",fill= col_fam)  +
    ylim(0,1800)+
    labs(title="(c)", x=NULL, y = "Count",
         subtitle = "News coverage by spider families")+
    theme_niwot()+
    theme(axis.text.x = element_text(angle = 70, hjust=1 ) )
)

# Panel d
col_gen <- c(rep("grey50",9),
            "darkorchid4",
            "grey50",
            "darkgoldenrod2",
            "darkorchid4",
            rep("brown3",2),"black")

# Loading symbol for venomous families
img <- png::readPNG("Danger_symbol.png")

(plot_genus <- ggplot(Bar_plot_genus, aes(x= Genus_2, y=N))+
    geom_bar(stat="identity",alpha=0.8,colour = "black", fill = col_gen)+
    ylim(0,1800)+
    labs(title="(d)", x=NULL, y = NULL,
         subtitle = "News coverage by spider genera\n[Danger symbols indicate genera with species of medical importance]")+
    
    annotation_custom(rasterGrob(img), xmin = 0.5, xmax = 1.5, ymin = 1600, ymax = 1800)   + 
    annotation_custom(rasterGrob(img), xmin = 9.5, xmax = 10.5, ymin = 1600, ymax = 1800)  +
    annotation_custom(rasterGrob(img), xmin = 12.5, xmax = 13.5, ymin = 1600, ymax = 1800) +
    annotation_custom(rasterGrob(img), xmin = 14.5, xmax = 15.5, ymin = 1600, ymax = 1800) +
    annotation_custom(rasterGrob(img), xmin = 15.5, xmax = 16.5, ymin = 1600, ymax = 1800) +
    # annotate("segment", 
    #          x    = c(1,10,13,15,16), 
    #          xend = c(1,10,13,15,16), 
    #          y    = 1800, 
    #          yend = 1600, 
    #          colour = "black", 
    #          size = 0.5, 
    #          arrow=arrow(ends = "last",angle = 45, length = unit(.2,"cm")))+
    # 
    theme_niwot()+
    theme(axis.text.x = element_text(angle = 70, hjust=1,
                                     face = c("plain",rep("italic",
                                                          nrow(Bar_plot_genus)-1))) )
)

#Arrange in a plate
plots <- list(map1, plot_lat, plot_family, plot_genus) ; grobs <- list() ; widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

pdf("Figure_3.pdf", width = 18, height = 10)

do.call("grid.arrange", c(grobs, nrow = 2, ncol = 2))

dev.off()

# Figure 4  ---------------------------------------------------------------

#re-ordering continents
db$Continent <- factor(db$Continent,
                       levels = c("Africa","Asia","Europe","N America","S America","Oceania"))

#Generating Bar plots

#Expert
bar_1 <- data.frame(table(db$TotalExpert,db$Continent))

(bar_p1 <-  ggplot(bar_1, aes(x=Var2,y=Freq, fill=Var1)) +
  
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("No", "Yes"),values=c("grey10", "turquoise3"))+
  labs(title="(a)", subtitle =  "Frequency of expert consultation (any)",x=NULL, y = "Frequency")+
  theme_niwot()+
  theme(legend.position = c(0.15, 0.8)))

#Spider expert
bar_2 <- data.frame(table(db$Expert_arachnologist,db$Continent))

(bar_p2 <-  ggplot(bar_2, aes(x=Var2,y=Freq, fill=Var1)) +
    
    geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
    geom_text(aes(label=Freq), vjust=-1, color="black",
              position = position_dodge(0.9), size=3.5)+
    
    scale_fill_manual("",labels=c("No", "Yes"),values=c("grey10", "turquoise3"))+
    labs(title="(b)", subtitle =  "Frequency of spider expert consultation",x=NULL, y = "Frequency")+
    theme_niwot()+
    theme(legend.position ="none"))

#Errors
bar_3 <- data.frame(table(ifelse(db$TotalError>0,1,0),db$Continent))

(bar_p3 <-  ggplot(bar_3, aes(x=Var2,y=Freq, fill=Var1)) +
    
    geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
    geom_text(aes(label=Freq), vjust=-1, color="black",
              position = position_dodge(0.9), size=3.5)+
    
    scale_fill_manual("",labels=c("No", "Yes"),values=c("grey10", "turquoise3"))+
    labs(title="(c)", subtitle =  "Presence of errors (any)",x=NULL, y = "Frequency")+
    theme_niwot()+
    theme(legend.position ="none"))

#Sensationalism
bar_4 <- data.frame(table(db$Sensationalism,db$Continent))

(bar_p4 <-  ggplot(bar_4, aes(x=Var2,y=Freq, fill=Var1)) +
  
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("No", "Yes"),values=c("grey10", "turquoise3"))+
  labs(title="(d)", subtitle =  "Sensationalistic content",x=NULL, y = "Frequency")+
  theme_niwot()+
  theme(legend.position = "none"))

#Arrange in a plate
plots <- list(bar_p1,bar_p2,bar_p3,bar_p4) ; grobs <- list() ; widths <- list()

for (i in 1:length(plots)){
  grobs[[i]] <- ggplotGrob(plots[[i]])
  widths[[i]] <- grobs[[i]]$widths[2:5]
}

maxwidth <- do.call(grid::unit.pmax, widths)

for (i in 1:length(grobs)){
  grobs[[i]]$widths[2:5] <- as.list(maxwidth)
}

#Arrange in a plate
pdf("Figure_4.pdf", width = 16, height = 12)

do.call("grid.arrange", c(grobs, nrow = 2, ncol = 2))

dev.off()

# Wordcloud (Figure 5) ----------------------------------------------------

#Extracting english news
db_english <- db[db$Lenguage == "English",]

#taking distinct titles
db_english_unique <- distinct(db_english, Title, .keep_all = TRUE) 

# A list of common and non-useful words, bundled with `tidytext`
data(stop_words)

# Searched words ("spider","spiders","bites", etc.) should also be omitted
searched_words <- c("false", # false black widow
                    "bananas", #banana spider
                    "bite",
                    "bites",
                    "bitten",
                    "spider",
                    "funnel", #funnel web spider
                    "latrodectus",
                    "loxosceles",
                    "spiders",
                    "sting",
                    "black",
                    "red",
                    "white",
                    "brown",
                    "back",
                    "redback",
                    "widow",
                    "widows",
                    "recluse",
                    "huntsman",
                    "tarantula",
                    "tarantulas")

# combine
stop_words <- data.frame(word=c(stop_words$word,searched_words))
                    
# remove all numbers from titles
db_english_unique$Title <- gsub('[0-9]+', '', db_english_unique$Title)

# extract the most common words in sensationalistic titles
Title_sensationalistic <- db_english_unique[db_english_unique$Sensationalism == 1,] %>%
  mutate(title = as.character(Title)) %>%
  unnest_tokens(output = title_word,
                input = Title) %>%
  anti_join(stop_words, by = c("title_word" = "word")) %>%
  count(title_word, sort = TRUE) 

#Setting colors
color_sens <- c("black", 
          "darkred",
          rep("black",6), #mum
          rep("darkorange",2),
          "black",
          "darkred", #horror
          "black",
          "darkred",
          "black",
          "darkred",
          "darkorange", #poisonous
          rep("black",7),
          "darkred", #terrifying
          rep("black",7),
          "darkred",#terrified
          "black",
          "darkred",
          "darkorange",
          "darkred",
          rep("black",2),
          "darkred",#hundreds
          rep("black",4),
          "darkred",
          "black",
          "darkorange",
          rep("black",13)
          )

dev.off()
(wordcloud1 <- Title_sensationalistic[1:60,] %>% with(wordcloud(words = title_word, 
                                         freq = n, 
                                         max.words = 60,
                                         scale=c(4,.2),
                                         random.color=FALSE, 
                                         ordered.colors = TRUE,
                                         colors = color_sens))
)


nrow(db_english_unique[db_english_unique$Sensationalism == 1,]) #sample size

# Extracting the most common words in non sensationalistic titles
Title_non <- db_english_unique[db_english_unique$Sensationalism == 0,] %>%
  mutate(title = as.character(Title)) %>%
  unnest_tokens(output = title_word,
                input = Title) %>%
  anti_join(stop_words, by = c("title_word" = "word")) %>%
  count(title_word, sort = TRUE) 

#Setting colors
color_non <- c(rep("black",11),
               "darkorange",#deadly
               rep("black",10),
               "darkred", #giant 
               rep("black",2),
               "darkorange", #kill
               rep("black",7),
               "darkred",#flesh
               rep("black",16),
               "darkorange",#kills
               rep("black",9)
               )

dev.off() 
(wordcloud2 <- Title_non[1:60,] %>% with(wordcloud(words = title_word, 
                                                  freq = n, 
                                                  max.words = 200,
                                                  scale=c(4,.2),
                                                  random.color=FALSE, 
                                                  ordered.colors = TRUE,
                                                  colors = color_non))
  
)

nrow(db_english_unique[db_english_unique$Sensationalism == 0,]) #sample size

## Note that Figure 5 has been assembled outside R, using Inkscape.

###############################################################

## Data validation:

###############################################################

# Percentager of reassessed news

#English
nrow(db[db$Quality_check == "yes" & db$Lenguage == "English",])/nrow(db[db$Lenguage == "English",])
#French
nrow(db[db$Quality_check == "yes" & db$Lenguage == "French",])/nrow(db[db$Lenguage == "French",])
#Italian
nrow(db[db$Quality_check == "yes" & db$Lenguage == "Italian",])/nrow(db[db$Lenguage == "Italian",])
#Spanish
nrow(db[db$Quality_check == "yes" & db$Lenguage == "Spanish",])/nrow(db[db$Lenguage == "Spanish",])

# Cohen's kappa analysis for calculating the correlation of two raters

# Loading the Database pre validation ----------------------------------------

db2 <- read.csv(file = "Data_spider_news_BeforeValidation.csv", sep = '\t', dec = '.', header = TRUE, as.is = FALSE)

db2 <- db2[db2$Lenguage == "English",] ; db2$Code <- droplevels(db2$Code)

# Extracting validated news
db1 <- db[db$Quality_check == "yes" & db$Lenguage == "English",] ; db1$Code <- droplevels(db1$Code)

# Match
db2 <- db2[db2$Code %in% db1$Code,]
db1 <- db1[db1$Code %in% db2$Code,]

# Sort
db1 <- db1 %>% arrange(factor(Code))
db2 <- db2 %>% arrange(factor(Code))

# Calculating inter-rater agreement ----------------------------------------

# Event
psych::cohen.kappa(x=cbind(db1$Bite,db2$Bite))
psych::cohen.kappa(x=cbind(db1$Death,db2$Death))

# Figure species
psych::cohen.kappa(x=cbind(db1$Figure_species,db2$Figure_species))

# Figure bite
psych::cohen.kappa(x=cbind(db1$Figure_bite,db2$Figure_bite))

# Experts
psych::cohen.kappa(x=cbind(db1$Expert_arachnologist,db2$Expert_arachnologist))
psych::cohen.kappa(x=cbind(db1$Expert_doctor,db2$Expert_doctor))
psych::cohen.kappa(x=cbind(db1$Expert_others,db2$Expert_others))

# Sensationalism
psych::cohen.kappa(x=cbind(db1$Sensationalism,db2$Sensationalism))

# Errors
psych::cohen.kappa(x=cbind(db1$Taxonomic_error,db2$Taxonomic_error))
psych::cohen.kappa(x=cbind(db1$Venom_error,db2$Venom_error))
psych::cohen.kappa(x=cbind(db1$Anatomy_error,db2$Anatomy_error))
psych::cohen.kappa(x=cbind(db1$Photo_error,db2$Photo_error))

#End of the script
