library(pophelper)
library(gridExtra)

###################
# PCA
###################

library("adegenet")
# read in from structure format

obj <- read.structure("~/Documents/GEOMAR/Teaching/Mar_pop_gen/pca/fheteroclitus_coastal.pops.stru",
                      n.ind=239,
                      n.loc= 15542,
                      onerowperind = TRUE,
                      col.lab = 1,
                      col.pop = 2,
                      row.marknames=1,
                      col.others=FALSE)

pop(obj) <- substr(indNames(obj), 1,2)
# replace the population names with the corresponding numbers from the map. Just to make it easier to interpret
pop(obj)<- gsub("TR", "1",pop(obj))
pop(obj)<- gsub("PC", "2",pop(obj))
pop(obj)<- gsub("HP", "3",pop(obj))
pop(obj)<- gsub("PP", "4",pop(obj))
pop(obj)<- gsub("PL", "5",pop(obj))
pop(obj)<- gsub("GA", "6",pop(obj))

x.af <- tab(obj, freq=TRUE, NA.method="mean")
pca.af <- dudi.pca(df = x.af, center = TRUE, # center by the mean
                   scale = FALSE, # don't scal/normalize,
                   scannf = FALSE, # don't show scree plot
                   nf = 5) # keep 5 pc's)

plot_data <- data.frame(individual = row.names(pca.af$li),
                        population = pop(obj),
                        PC1 = pca.af$li[,1],
                        PC2 = pca.af$li[,2]
)
plot_data$population <- factor(plot_data$population, levels = c("1","2","3","4","5","6"))


ggplot(plot_data, aes(x=PC1, y=PC2, fill=population)) +
  geom_point(size=3, shape=21) +
  theme_classic()

eig.perc <- 100*pca.af$eig/sum(pca.af$eig)
barplot(eig.perc[1:20])

###############
# structure
###############

files <- list.files(path="~/Documents/GEOMAR/Teaching/Mar_pop_gen/structure",
                    full.names=TRUE)
slist <- readQ(files=files, filetype="structure", indlabfromfile=T)
head(slist[[15]])

names(attributes(slist[[2]]))
qlist <- (tabulateQ(slist))
sr1 <- summariseQ(qlist)
#evannoMethodStructure(data=sr1)

p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)

slist <- alignK(slist)
length(slist)
mergedS <- (mergeQ(slist))

# make group labels:
labs <- data.frame(location = substring(rownames(slist[[15]]),1,2))

p1 <- plotQ(mergedS,returnplot=T,exportplot=F,basesize=11,
            sortind="Cluster1", grplab=labs)
grid.arrange(p1$plot[[1]], p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],p1$plot[[7]],p1$plot[[7]])

