######################################################################
# re-analysis of Freeman et al. 2022 Science paper
# by Dingliang Xing xingdingliang@gmail.com, July 29, 2022
######################################################################
library(raster)
library(sf)
library(lavaan)
library(ggplot2)
library(dplyr)
library(glue)
library(here)
library(arrow)
theme_set(theme_classic())
# part 1: data ------------------------------------------------------
## data preparation -------------------------------------------------

# this section is mostly adapted from Freeman's original paper.

# set working directory to folder that contains files with species coded as forest or not.
# this is the path in the unzipped file that was downloaded from Freeman's paper.
setwd("Part2_code_data/Part2/files/") 

## 
mydata <- read.csv("elevation_dataframe.csv")
mydata <- mydata[mydata$include == "yes",] # subset to species that live in forest
dim(mydata) # 5397 rows
length(unique(mydata$common_name)) # 2879 unique species


# calculate elevational range size for each species - using the frequency method to define elevational ranges
mydata$breadth <- mydata$upper_frequency - mydata$lower_frequency

# Xing added: calculate elevational position for each species
mydata$position <- (mydata$upper_frequency + mydata$lower_frequency)/2


# make a dataframe where each region is a different row, with several attributes summarized for each region (e.g., species richness, elevational range size, etc)

# Xing added rposition

summarydata <- mydata %>%   # Specify data frame
    group_by(region) %>%      # Specify group indicator
    dplyr::summarise(sprichness = n(), breadth = mean (breadth), elev_min = min(lower_frequency), elev_max = max(upper_frequency), sampling = median(n_pos_obs), sampling_mean = mean(n_pos_obs), expanse=elev_max-elev_min, position_median = median(position))

# add latitudes for regions from a different csv; bio4 in this csv
latitudes <- read.csv("spatial_info_for_regions.csv")

latitudes$region <- tolower(latitudes$region) # so that the two dataframes match case

summarydata2 <- left_join(summarydata, latitudes, by = "region")

# add climate change velocity for regions from a different csv
velocities <- read.csv("spatial_info_for_regions_velocity.csv")
velocities$region <- tolower(velocities$region) # so that the two dataframes match case

summarydata2 <- left_join(summarydata2, velocities, by = "region")

# add sampling completeness metric for each region based on sampling completeness analyses, again from a different csv
completeness_df <- read.csv("completeness_df.csv")
summarydata2 <- left_join(summarydata2, completeness_df, by = "region")


# Xing added: calculate and add area for each region
summarydata2$area = NA
regionlist = dir('shp_files_for_MSM/', '*.shp')
for(i in regionlist[-7]){
    ind_i = which(summarydata2$region==tolower(gsub('-polygon.shp','',i)))
    region_i = sf::st_read(paste0('shp_files_for_MSM/',i))
    area_i = sf::st_area(region_i)/1e6 # area in km2
    summarydata2$area[ind_i] = area_i
}
# the spherical geometry of the polygon for Cascades is invalid. estimate the area by switching spherical geometry off:
sf::sf_use_s2(FALSE)
i=regionlist[7]
ind_i = which(summarydata2$region==tolower(gsub('-polygon.shp','',i)))
region_i = sf::st_read(paste0('shp_files_for_MSM/',i))
area_i = sf::st_area(region_i)/1e6 # area in km2
summarydata2$area[ind_i] = area_i
sf::sf_use_s2(TRUE)
summarydata2$rposition_median = (summarydata2$position_median -  summarydata2$elev_min)/summarydata2$expanse

write.csv(summarydata2, 'eBirdSummarydata.csv', row.names = F)

## Fig 1 ------------------------------------------------------------
### Fig 1A -------
# get elevation data for the three countries in this region - Ecuador Colombia & Peru 
# getData is from the raster package
alt = getData("alt", country = "ECU", path = tempdir())
alt2 = getData("alt", country = "COL", path = tempdir())
alt3 = getData("alt", country = "PER", path = tempdir())

# put these together
elev <- merge(alt, alt2, alt3)
plot(elev) # ok it worked

#convert the raster to points for plotting
elev.p <- rasterToPoints(elev)

#Make the points a dataframe so that ggplot can use it
df <- data.frame(elev.p)

#Make appropriate column headings for the dataframe
colnames(df) <- c("Longitude", "Latitude", "Elev")


#make a bounding box by subsetting this dataframe

df_plotting <- df[df$Longitude > -80.5,]
df_plotting <- df_plotting[df_plotting$Longitude < -72.5,]
df_plotting <- df_plotting[df_plotting$Latitude < 5,]
df_plotting <- df_plotting[df_plotting$Latitude > -5,]

# Xing added: load and add the shapfiles of the two regions on the map
ChocoShp = sf::st_read('shp_files_for_MSM/Choco-polygon.shp')
AmazonShp = sf::st_read('shp_files_for_MSM/C_Tropical_Andes-polygon.shp')

elev_Choco = data.frame(rasterToPoints(mask(crop(elev,ChocoShp), ChocoShp)))
elev_Amazon = data.frame(rasterToPoints(mask(crop(elev,AmazonShp), AmazonShp)))

# load the sampling locations
checklistsc <- glue::glue("choco_checklists.parquet") %>%
    here::here("Part1_code_data/Code_MSM/Part1/data", "ebird", .) %>%
    arrow::read_parquet() %>%
    filter(!is.na(elevation), elevation >= 0)
checklistsa <- glue::glue("c_tropical_andes_checklists.parquet") %>%
    here::here("Part1_code_data/Code_MSM/Part1/data", "ebird", .) %>%
    arrow::read_parquet() %>%
    filter(!is.na(elevation), elevation >= 0)


ggplot()+
    geom_raster(data=df_plotting, 
                aes(y=Latitude, x=Longitude,fill=sqrt(Elev))) +
    scale_fill_gradient("Elev (m)",
                        breaks=c(0,20,40,60),
                        labels=c(0,400,1600,3600),
                        low = "gray80", high = "gray28") +
    coord_equal()+
    # Xing: add the sampling locations and polygons for the two regions
    geom_point(data=checklistsc, aes(x=longitude, y=latitude), size=0.1)+
    geom_point(data=checklistsa, aes(x=longitude, y=latitude), size=0.1)+
    geom_sf(data=ChocoShp,fill=rgb(0,0,0,.1))+
    geom_sf(data=AmazonShp,fill=rgb(0,0,0,.1))+
    annotate(geom="text", x=-76, y=-4.5, label="Amazon",
             size=6)+
    annotate(geom="text", x=-79.2, y=3.5, label="Chocó",
             size=6)+
    theme(legend.position=c(.95,.95),
          legend.justification = c("right", "top"))

ggsave(here("fig1a.svg"), device = 'svg', height = 4, width = 3.5)

### Fig 1B -----
mydata2 = left_join(mydata, summarydata2[,c('region', 'elev_min', 'expanse')])

comparison <- mydata2[mydata2$region == "choco" | mydata2$region == "c_tropical_andes",]
comparison$region = as.factor(comparison$region)
comparison$region = relevel(comparison$region, 'choco')

areachoco = st_area(ChocoShp)/1e6/1e3 # unit: km2*1000
areaamazon = st_area(AmazonShp)/1e6/1e3 # unit: km2*1000

area_comp = data.frame(region=c('Chocó', 'Amazon'),
                       area=c(areachoco, areaamazon)) 
area_comp$area = as.numeric(area_comp$area)
area_comp$region = as.factor(area_comp$region)
area_comp$region = relevel(area_comp$region, 'Chocó')

ggplot(data = area_comp, aes(x=region, y=area, fill=region))+
    geom_bar(stat = "identity") +
    scale_y_continuous(expression("Area (1000 "~km^2~")")) +
    scale_x_discrete("")+
    theme(legend.position = "none")

ggsave(here("fig1b.svg"), device = 'svg', height = 4, width = 3)

### Fig 1C -------
comparison$rposition=(comparison$position-comparison$elev_min)/comparison$expanse

t.test(comparison$rposition~comparison$region)
# 
# Welch Two Sample t-test
# 
# data:  comparison$rposition by comparison$region
# t = 3.5746, df = 1060.1, p-value = 0.0003666
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#     0.02087672 0.07168838
# sample estimates:
#     mean in group choco mean in group c_tropical_andes 
# 0.3547266                      0.3084440 
effsize::cohen.d(comparison$rposition[comparison$region=='choco'],
                 comparison$rposition[comparison$region=='c_tropical_andes'])
# 
# Cohen's d
# 
# d estimate: 0.2154626 (small)
# 95 percent confidence interval:
#      lower      upper 
# 0.09697163 0.33395364 

comparison %>% 
    ggplot(aes(x=region,y=rposition, fill=region)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(position=position_jitterdodge(jitter.width = 1, jitter.height = .001, dodge.width = 0),alpha=0.1)+
    scale_y_continuous("Relative elevational position") +
    scale_x_discrete("", labels = c("Chocó","Amazon"))+
    theme(legend.position = 'none')

ggsave(here("fig1c.svg"), device = 'svg', height = 4, width = 3)



# part 2: multiple regression ---------------------------------------
## Freeman's regression model--------- 
lm1 <- lm(breadth ~ completeness_median + expanse + bio4 + sprichness + cc_velocity, data = summarydata2)
summary(lm1)

# excluding the four most species-rich regions (all from Andes and have low median species elevational position) will give completely different results;
# the effect of sprichness is no longer significant. bio4 is also not significant, though.
lm1sub <- lm(breadth ~ completeness_median + expanse + bio4 + sprichness + cc_velocity, data = subset(summarydata2, sprichness<400))
summary(lm1sub)
# sprichness is correlated with position:
with(summarydata2, plot(rposition_median, sprichness))
with(summarydata2, cor.test(rposition_median, sprichness))
# Pearson's product-moment correlation
# 
# data:  rposition_median and sprichness
# t = -5.724, df = 29, p-value = 3.405e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.8605776 -0.5040870
# sample estimates:
#        cor 
# -0.7283335 

## our new regression model---------
# sprichness and bio4 are highly right skewed:
with(summarydata2, hist(sprichness))
with(summarydata2, hist(bio4))
# take log on both variables:
summarydata2$sprichness_log=log(summarydata2$sprichness)
summarydata2$bio4_log=log(summarydata2$bio4)

# our multiple regression model: 
lm1_my <- lm(breadth ~ completeness_median + expanse + bio4_log + sprichness_log + cc_velocity + rposition_median, data = summarydata2)

# there is no multicollinearity problem:
car::vif(lm1_my)
# completeness_median             expanse            bio4_log 
# 1.562052            2.829639            1.931563 
# sprichness_log         cc_velocity    rposition_median 
# 3.625443            1.414690            1.865169 

# comparison between the new model and their original model:
# my model has larger adjusted R2
summary(lm1)
# Adjusted R-squared:  0.686 
summary(lm1_my)
# Adjusted R-squared:  0.7204

# my model has lower AIC
AIC(lm1)
# 419.18
AIC(lm1_my)
# 416.3098

# the residuals of my model are more normally distributed
shapiro.test(lm1$residuals)
# 
# Shapiro-Wilk normality test
# 
# data:  lm1$residuals
# W = 0.93007, p-value = 0.04404
shapiro.test(lm1_my$residuals)
# 
# Shapiro-Wilk normality test
# 
# data:  lm1_my$residuals
# W = 0.96845, p-value = 0.4775

## Fig 2B and 2C ---------
sensemakr::partial_f2(lm1_my, covariates = "bio4_log")
# 0.2343108    
sensemakr::partial_f2(lm1_my, covariates = "sprichness_log")
# 0.02221261    

rpos=summarydata2$rposition_median

visreg::visreg(lm1_my, xvar = "sprichness_log", gg = T) +
    geom_point(aes(colour=rpos))+
    scale_colour_gradientn(colours = terrain.colors(30)[1:20])+
    guides(colour='none',size='none')+
    scale_x_continuous("Regional species richness (log)") +
    scale_y_continuous("Partial residual elevational range size (m)", limits = c(0, 2000)) 
ggsave(here("fig2b.svg"), device = 'svg', height = 3.5, width = 3.5)

visreg::visreg(lm1, xvar = "sprichness", xtrans=log, gg = T) +
    geom_point(aes(colour=rpos))+
    scale_colour_gradientn(colours = terrain.colors(30)[1:20])+
    guides(colour='none',size='none')+
    scale_x_continuous("") +
    scale_y_continuous("", limits = c(0, 2000))+
    theme(plot.background = element_rect(fill = "grey"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
ggsave(here("fig2b_insert.svg"), device = 'svg', height = 1.7, width = 1.7)

visreg::visreg(lm1_my, xvar = "bio4_log", gg = T) +
    geom_point(aes(colour=rpos))+
    scale_colour_gradientn(name='Median species\nelevational position',colours = terrain.colors(30)[1:20])+
    guides(colour='none')+
    scale_x_continuous("Temperature seasonality (log)") +
    scale_y_continuous("Partial residual elevational range size (m)", limits = c(0, 2000))+
    theme(legend.position=c(.95,.02),
          legend.title=element_text(size=7),
          legend.text=element_text(size=7),
          legend.direction='horizontal',
          legend.justification = c("right", "bottom"))+
    guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))
ggsave(here("fig2c.svg"), device = 'svg', height = 3.5, width = 3.5)

## our results are robust------

# taking log on sprichness alone gives qualitatively similar results
lm1_logRich <- lm(breadth ~ completeness_median + expanse + bio4 + sprichness_log + cc_velocity, data = summarydata2)
summary(lm1_logRich)

# taking log on both sprichness and bio4 also gives qualitatively similar results
lm1_logRich_logBio4 <- lm(breadth ~ completeness_median + expanse + bio4_log + sprichness_log + cc_velocity, data = summarydata2)
summary(lm1_logRich_logBio4)

# adding the position effect gives my final model that has a much better model-data fit
lm1_my <- lm(breadth ~ completeness_median + expanse + bio4_log + sprichness_log + cc_velocity + rposition_median, data = summarydata2)
summary(lm1_my)

# further taking log on expanse does not change the results qualitatively
lm1_my2 <- lm(breadth ~ completeness_median + log(expanse) + bio4_log + sprichness_log + cc_velocity + rposition_median, data = summarydata2)
summary(lm1_my2)

# part 3: path analysis ---------------------------------------------

# use scaled variables for easier interpretation of coefficients

summarydata2$breadth_scaled <- scale(summarydata2$breadth)
summarydata2$expanse_scaled <- scale(summarydata2$expanse)
summarydata2$sprichness_scaled <- scale(summarydata2$sprichness)
summarydata2$sprichness_log_scaled <- scale(summarydata2$sprichness_log)
summarydata2$bio4_scaled <- scale(summarydata2$bio4)
summarydata2$bio4_log_scaled <- scale(summarydata2$bio4_log)
summarydata2$cc_velocity_scaled <- scale(summarydata2$cc_velocity)
summarydata2$completeness_median_scaled <- scale(summarydata2$completeness_median)
summarydata2$rposition_median_scaled = scale(summarydata2$rposition_median)
summarydata2$area_log_scaled=scale(log(summarydata2$area))



# Freeman et al.'s model:
model <- "breadth_scaled ~ expanse_scaled + sprichness_scaled + bio4_scaled + completeness_median_scaled
          sprichness_scaled ~  expanse_scaled + bio4_scaled + cc_velocity_scaled"
# my new model:
model_my <- "breadth_scaled ~ expanse_scaled + r2y*sprichness_log_scaled + s2y*bio4_log_scaled + completeness_median_scaled + cc_velocity_scaled + rposition_median_scaled
          sprichness_log_scaled ~ expanse_scaled + s2r*bio4_log_scaled + cc_velocity_scaled + area_log_scaled
          s2r2y:=s2r*r2y # indirect effect of bio4 on breadth
          s2ytot:=s2y+s2r2y # total effect of bio4 on breadth
"
fit= cfa(model, data=summarydata2)#Freeman et al.'s
fit_my= cfa(model_my, data=summarydata2)# My
summary(fit, fit.measures = TRUE, standardized=T, rsquare=T)
summary(fit_my, fit.measures = TRUE, standardized=T, rsquare=T)
# lavaan 0.6-7 ended normally after 29 iterations
# 
# Estimator                                         ML
# Optimization method                           NLMINB
# Number of free parameters                         12
# 
# Number of observations                            31
# 
# Model Test User Model:
#     
#     Test statistic                                 0.991
# Degrees of freedom                                 3
# P-value (Chi-square)                           0.803
# 
# Model Test Baseline Model:
#     
#     Test statistic                                90.063
# Degrees of freedom                                13
# P-value                                        0.000
# 
# User Model versus Baseline Model:
#     
#     Comparative Fit Index (CFI)                    1.000
# Tucker-Lewis Index (TLI)                       1.113
# 
# Loglikelihood and Information Criteria:
#     
#     Loglikelihood user model (H0)                -42.422
# Loglikelihood unrestricted model (H1)        -41.926
# 
# Akaike (AIC)                                 108.843
# Bayesian (BIC)                               126.051
# Sample-size adjusted Bayesian (BIC)           88.665
# 
# Root Mean Square Error of Approximation:
#     
#     RMSEA                                          0.000
# 90 Percent confidence interval - lower         0.000
# 90 Percent confidence interval - upper         0.189
# P-value RMSEA <= 0.05                          0.821
# 
# Standardized Root Mean Square Residual:
#     
#     SRMR                                           0.013
# 
# Parameter Estimates:
#     
#     Standard errors                             Standard
# Information                                 Expected
# Information saturated (h1) model          Structured
# 
# Regressions:
#     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# breadth_scaled ~                                                             
# expns_sc                 1.006    0.146    6.901    0.000    1.006    1.011
# sprchn__ (r2y)          -0.134    0.156   -0.858    0.391   -0.134   -0.135
# b4_lg_sc (s2y)           0.318    0.120    2.645    0.008    0.318    0.320
# cmpltn__                 0.162    0.105    1.538    0.124    0.162    0.163
# cc_vlct_                -0.203    0.102   -1.995    0.046   -0.203   -0.204
# rpstn_m_                 0.466    0.114    4.099    0.000    0.466    0.468
# sprichness_log_scaled ~                                                      
# expns_sc                 0.667    0.097    6.870    0.000    0.667    0.667
# b4_lg_sc (s2r)          -0.372    0.103   -3.601    0.000   -0.372   -0.372
# cc_vlct_                 0.181    0.098    1.849    0.064    0.181    0.181
# ar_lg_sc                 0.229    0.092    2.484    0.013    0.229    0.229
# 
# Variances:
#     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# .breadth_scaled    0.216    0.055    3.937    0.000    0.216    0.226
# .sprchnss_lg_sc    0.245    0.062    3.937    0.000    0.245    0.253
# 
# R-Square:
#     Estimate
# breadth_scaled    0.774
# sprchnss_lg_sc    0.747
# 
# Defined Parameters:
#     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# s2r2y             0.050    0.060    0.835    0.404    0.050    0.050
# s2ytot            0.368    0.114    3.220    0.001    0.368    0.370