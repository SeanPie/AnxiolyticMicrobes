library(readxl)
library(ggplot2)
library(eulerr)
library(dplyr)
library(tidyr)

setwd("/Users/sean/Desktop/Fall 21/BIOL 806/Lab")

datadir <- "data/Final"

colors <- c("Control" = "#f46d43", "Low" = "#66c2a5" , "High" = "#5e4fa2" ) # orange, aqua, dark-blue

### Figures for Elevated Plus Maze

ep_post <- read_excel(path= paste0(datadir,"/","Elev-plus-post-exposure.xlsx"),
                      sheet=1, range="A1:AA55", col_names = TRUE)
ep_post <- as.data.frame(ep_post)


ep_post$Treatment <- factor(ep_post$Treatment, levels = c("Control", "Low" , "High"), ordered = TRUE )

ep_post$non_Open_arm_time <- ep_post$Closed_arm_time + ep_post$Centre_time

ep_post$time_center_and_open_arm <- ep_post$Centre_time + ep_post$Open_arm_time



a <- ggplot(data = ep_post, aes(color = Treatment, x = Treatment, y = Open_arm_time)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  facet_grid(. ~ Sex) +
  labs(title = "Elevated Plus Maze Time Spent in Anxienty Inducing Areas", x = "Treatment", y = "Time in Open Arms (s)") +
  theme_classic()


b <- ggplot(data = ep_post, aes(color = Treatment, x = Treatment, y = time_center_and_open_arm)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  facet_grid(. ~ Sex) +
  labs(title = "Elevated Plus Maze Time Spent Including Center", x = "Treatment", y = "Time in Open Arms and Center (s)") +
  theme_classic()

### Figures for Open Field Test

of_post <- read_excel(path= paste0(datadir,"/","Open-field-post-exposure.xlsx"),
                      sheet=1, range="A1:AE55", col_names = TRUE)

of_post <- as.data.frame(of_post)

of_post$Treatment <- factor(of_post$Treatment, levels = c("Control", "Low" , "High"), ordered = TRUE )


c <- ggplot(data = of_post, aes(color = Treatment, x = Treatment, y = Centre_time)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  facet_grid( ~ Sex) +
  theme_classic() +
  labs(title = "Open Field Test for Time Spent in Anxiety Inducing Areas", x = "Treatment", y = "Time in Field Center (s)")

#Filtering for female mice only in Open Field Test
sel <- which(of_post$Sex == "female")
of_fem <- of_post[sel, c("ID", "Treatment", "Centre_time")]

d <- ggplot(data = of_fem, aes(color = Treatment, x = Treatment, y = Centre_time)) +
  geom_boxplot() +
  scale_color_manual(values = colors) +
  theme_classic() +
  labs(title = "Open Field Test for Only Females in Center", x = "Treatment", y = "Time in Field Center (s)")

a
b
c
d

raw_otu <- read_excel(path= paste0(datadir,"/","absoluteabundance.xlsx"),
                      sheet=1, range="A1:LN10485", col_names = TRUE)

raw_otu <- as.data.frame(raw_otu)

#Remove "-16S-V3-V4" from column names and change the name of the first column
names(raw_otu) <- gsub(pattern = "-16S-V3-V4",x = names(raw_otu), replacement = "")
names(raw_otu)[1] <- "OTU_ID"

#Remove OTUs that weren't bacteria taxa
raw_otu <- raw_otu[ grep("k__Archaea", raw_otu$`Consensus Lineage`, invert = TRUE) , ]
raw_otu <- raw_otu[ grep("Unclassified", raw_otu$`Consensus Lineage`, invert = TRUE) , ]

#Remove OTUs for chloroplast, streptophyta, and mitochondria
raw_otu <- raw_otu[ grep("c__Chloroplast", raw_otu$`Consensus Lineage`, invert = TRUE) , ]
raw_otu <- raw_otu[ grep("o__Streptophyta", raw_otu$`Consensus Lineage`, invert = TRUE) , ]
raw_otu <- raw_otu[ grep("f__mitochondria", raw_otu$`Consensus Lineage`, invert = TRUE) , ]

#Select only columns for mice fecal samples at week 1 (T01) and week 17 (T16)
clean_otu <- raw_otu %>% select(matches("m")) %>% select(matches("T"))
clean_otu$OTU_ID <- raw_otu$OTU_ID 
clean_otu <- clean_otu %>% select(OTU_ID, everything())

#Filters OTUs that didn't have atleast 50 sequence reads in more than one mice fecal sample.
clean_otu <- clean_otu[rowSums(clean_otu > 50) > 1, ]


### Figures for total Number of OTUs Pre and Post Dust Exposures

#Filtering and created a table for the total number of OTUs in mice pre dust exposure
otu_pre_exposure <- clean_otu %>% select(matches("T01"))
otu_pre_exposure <- otu_pre_exposure[, sort(names(otu_pre_exposure))]
pre_otu <- as.data.frame(colSums(otu_pre_exposure != 0))
pre_otu <- cbind(rownames(pre_otu), data.frame(pre_otu, row.names=NULL))
names(pre_otu) <- c("Mouse_ID", "Pre_otu")
pre_otu[1] <- gsub(pattern = "T01",x = (pre_otu$Mouse_ID), replacement = "")
ep_post <- ep_post %>% arrange(ID)
pre_otu$Treatment <- ep_post$Treatment
pre_otu$Sex <- ep_post$Sex


#Filtering and created a table for the total number of OTUs in mice post dust exposure
otu_post_exposure <- clean_otu %>% select(matches("T16"))
otu_post_exposure <- otu_post_exposure[, order(names(otu_post_exposure))]
post_otu <- as.data.frame(colSums(otu_post_exposure != 0))
post_otu <- cbind(rownames(post_otu), data.frame(post_otu, row.names=NULL))
names(post_otu) <- c("Mouse_ID", "Post_otu")
post_otu[1] <- gsub(pattern = "T16",x = (post_otu$Mouse_ID), replacement = "")
ep_post <- ep_post %>% arrange(ID)
post_otu$Treatment <- ep_post$Treatment
post_otu$Sex <- ep_post$Sex

df1 <- aggregate(pre_otu$Pre_otu, list(pre_otu$Treatment), FUN=mean)
names(df1) <- c("Treatment", "Pre")
df2 <- aggregate(post_otu$Post_otu, list(post_otu$Treatment), FUN=mean)
names(df2) <-  c("Treatment", "Post")
df3 <- full_join(df1, df2)
df3 <- gather(df3, key="Exposure", value="OTUs", 2:3)

e <- ggplot(data = df3, aes(color = Treatment, group = Treatment, x = Exposure, y = OTUs)) +
  scale_x_discrete(limits = rev) +
  scale_color_manual(values = colors) +
  geom_point() +
  geom_line() +
  theme_classic()

e

#Filtering for pre and post exposure for only female mice
sel <- which(pre_otu$Sex == "female")
pre_fem <- pre_otu[sel, c("Treatment", "Pre_otu")]
pre_fem <- aggregate(pre_fem$Pre_otu, list(pre_fem$Treatment), FUN=mean)
names(pre_fem) <- c("Treatment", "Pre")

sel <- which(post_otut$Sex == "female")
post_fem <- post_otu[sel, c("Treatment", "Post_otu")]
post_fem <- aggregate(post_fem$Post_otu, list(post_fem$Treatment), FUN=mean)
names(post_fem) <- c("Treatment", "Post")
df4 <- full_join(pre_fem, post_fem)
df4 <- gather(df4, key="Exposure", value="OTUs", 2:3)


f <- ggplot(data = df4, aes(color = Treatment, group = Treatment, x = Exposure, y = OTUs)) +
  scale_x_discrete(limits = rev) +
  scale_color_manual(values = colors) +
  geom_point() +
  geom_line() +
  theme_classic()

f








