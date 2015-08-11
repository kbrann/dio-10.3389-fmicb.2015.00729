require(grid)
library(phyloseq)
require(ggplot2)
require(vegan)
require(RColorBrewer)
require(gridExtra)
set.seed(123)

#get our water data
#import qiime biom with taxonomy and tre file
taxvec1 = c("k__Bacteria", "p__Firmicutes", "c__Bacilli", "o__Bacillales", "f__Staphylococcaceae")
mydata<-import_biom("otu_table_mc2_w_tax.biom","rep_set.tre",parseFunction=parse_taxonomy_greengenes)
parse_taxonomy_greengenes(taxvec1)
#import mapping file
qiimedata <- import_qiime_sample_data("mapping.txt")
#merge them together
data <- merge_phyloseq(mydata, qiimedata)
#Finally, the following is the remaining set of preprocessing steps that was applied to the GlobalPatterns OTU counts prior to creating the figures in the main phyloseq manuscript.
#Remove taxa not seen more than 3 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
data<-prune_samples(c("CCR.13d2",
                      "CCR.13d1",
                      "CCR.13N1",
                      "CCR.13N2",
                      "CCR.13O.2A",
                      "CCR.13O.2B",
                      "CCR.13O.1A",
                      "CCR.13O.1B",
                      "CCR.13O.1C",
                      "CCR.13S.3A",
                      "CCR.13S.3B",
                      "CCR.13S.3C",
                      "CCR.13S.2A",
                      "CCR.13S.2B",
                      "CCR.13S.2C",
                      "CCR.13A3",
                      "CCR.13A2",
                      "CCR.13A1",
                      "CCR13D3",
                      "CCR13D2",
                      "CCR13D1",
                      "CCR13N3",
                      "CCR13N2",
                      "CCR13N1",
                      "CCR.13O3",
                      "CCR.13O2",
                      "CCR.13O1",
                      "CCR.13SF3",
                      "CCR.13S3",
                      "CCR13SL3",
                      "CCRS.13A3",
                      "CCRS.13A2",
                      "CCRS.13A1",
                      "CCRS.13d2",
                      "CCRS.13d1",
                      "CCRS.13N3",
                      "CCRS.13N2",
                      "CCRS.13N1",
                      "CCRS.13O3",
                      "CCRS.13O2",
                      "CCRS.13O1",
                      "CCRS.13S3",
                      "CCRS.13S2","CCR7","CCR8","CCR9","CCR10","CCR11"
),data)

data_filter= filter_taxa(data, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
ntaxa(data)
ntaxa(data_filter)
#Standardize abundances to the median sequencing depth

total = median(sample_sums(data_filter))
standf = function(x, t=total) round(t * (x / sum(x)))
data_trans = transform_sample_counts(data_filter, standf)

#Filter the taxa using a cutoff of 2.0 for the Coefficient of Variation

data_trans_cv = filter_taxa(data_trans, function(x) sd(x)/mean(x) > 2.0, TRUE)
ntaxa(data_trans_cv)

#Ordinate NMDS using Bray, use stress plot to find the best dissimilarit/distance
nmds.bray <- ordinate(data,"NMDS","bray") 
stressplot(nmds.bray) 
as.data.frame(nmds.bray$points)

environ<-data.frame(sample_data(data)[,2:16])
str(environ)
##site data
sites <- data.frame(scores(nmds.bray, choices=c(1,2),display = c("sites"))) #dataframe of species scoes for plotting
head(sites)
sites$Location <- sample_data(data)$Location #otherwise factor doesn't drop unused levels and it will throw an error
sites$Type <- sample_data(data)$Type
sites$Season <- sample_data(data)$Season

#envfit
nmds.bray.envfit <- envfit(as.data.frame(nmds.bray$points), 
                           env = environ,na.rm=TRUE, perm = 999) #standard envfit
plot(nmds.bray)
plot(nmds.bray.envfit,p=0.05)

#ellipses
# function for ellipsess - just run this, is used later
#taken from the excellent stackoverflow Q+A: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the management factor
df_ell.dune.management <- data.frame() #sets up a data frame before running the function.
for(g in levels(sites$Type)){
  df_ell.dune.management <- rbind(df_ell.dune.management, cbind(as.data.frame(with(sites [sites$Type==g,],
                                                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                                                ,Type=g))
}

# data for labelling the ellipse
NMDS.mean.biotrap=aggregate(sites[ ,c("NMDS1", "NMDS2")],
                            list(group = sites$Type), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(sites[,c("NMDS1", "NMDS2")],
                    list(group = sites$Type), mean)


## finally plotting.
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
#display.brewer.pal(8,"Dark2")

nmds2 <- ggplot(data = sites, aes(x = NMDS1, y = NMDS2))+ #sets up the plot. brackets around the entire thing to make it draw automatically
  geom_point(aes(x = NMDS1, y = NMDS2, shape=Type),size = 6) + #puts the site points in from the ordination, shape determined by site, size refers to size of point
  geom_path(data = df_ell.dune.management, aes(x = NMDS1, y = NMDS2, group = Type)) + #this is the ellipse, seperate ones by type. If you didn't change the "alpha" (the shade) then you need to keep the "group  annotate("text",x = -0.70,y = 0.5,label="Filter",size=6) + #labels for the centroids - I haven't used this since we have a legend. but you could also dithc the legend, but plot will get v messy
  theme_bw(base_size = 15) +
  scale_manual_color() +
  theme(axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(), axis.title.y=element_blank()
  )
            

ggsave(file="NMDS_plot.tiff",plot=nmds2,width=10,height=7,units="in",dpi=300)



