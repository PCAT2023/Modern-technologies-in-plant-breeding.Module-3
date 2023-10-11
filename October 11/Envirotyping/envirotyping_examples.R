################ Install EnvRtype ################

install.packages("devtools")
devtools::install_github("allogamous/EnvRtype", force = TRUE)
require("EnvRtype")

setwd("I:/My Drive/Skoltech OneDrive/Courses/MPBW_23")

################ Practical use of get_weather ################

env.data = get_weather(env.id = 'Ust-Labinsk', country = 'RUS',
                       lat = 45.262088,lon = 39.678406,
                       start.day = '2000-03-01',end.day = '2021-04-01')

head(env.data)


################ Practical use of SummaryWTH ################ 
summaryWTH(env.data = env.data, env.id = 'env', days.id = 'daysFromStart',statistic = 'mean')


################ Basic use of env_typing for typing temperature from 2000 to 2020 ################ 

env.data = get_weather(env.id = 'Ust-Labinsk',country = 'RUS',
                       lat = 45.262088,lon = 39.678406,variables.names = 'T2M',
                       start.day = '2000-03-01',end.day = '2020-03-01')

card = list(T2M=c(0,8,15,28,40,45,Inf)) # a list of vectors containing empirical and cardinal thresholds
env_typing(env.data = env.data,env.id = 'env', var.id = 'T2M', cardinals = card)

################ Basic use of env_typing for more than one variable ################
var = c("PRECTOT", "T2MDEW") # variables
env.data = get_weather(env.id = 'Ust-Labinsk',country = 'RUS',
                       lat = 45.262088,lon = 39.678406,variables.names = var,
                       start.day = '2000-03-01',end.day = '2020-03-01')
card = list(PRECTOT = c(0,5,10,25,40,100), T2MDEW = NULL) # cardinals and data-driven limits
env_typing(env.data = env.data,env.id = 'env', var.id = var, cardinals = card)

################ Basic use of env_typing for more than one variable ################
data("maizeWTH") # toy set of environmental data
var = c("PRECTOT", "T2MDEW", "T2M_MAX", "T2M_MIN") # variables
W = W_matrix(env.data = maizeWTH[maizeWTH$daysFromStart < 100,],
             var.id=var, statistic="mean", by.interval=TRUE)
dim(W)


################ Remote Sensing for Several Places ################
env = c('GOI','TEX','BRI','MON','LOS','PON','CAL','PAL','DAV')
lat = c(-16.67,19.25,-27.47,43.61,14.170,6.294,3.261,-10.168,38.321)
lon = c(-49.25,-99.50,153.02,3.87,121.241,2.361,-76.312,-48.331,-121.442)
#cou = c('BRA', 'MEX', 'AUS', 'FRA', 'PHI', 'BEN', 'COL', 'BRA', 'USA')
start = c('2020-03-15','2019-05-15','2018-09-15',
          '2017-06-18','2017-05-18','2016-07-18',
          '2017-11-18','2017-12-18','2018-07-18')
end = c('2020-04-15','2019-06-15','2018-10-15',
        '2017-07-18','2017-06-18','2016-08-18',
        '2017-12-18','2018-01-18','2018-08-18')
env.data = get_weather(env.id = env, lat = lat, lon = lon, start.day = start, end.day = end)


################ Discovering ETs and similarity among locations ################

ET = env_typing(env.data = env.data,env.id = 'env',var.id = var,format = 'wide')
#ET = env_typing(env.data = env.data,env.id = 'env',var.id = 'T2M',format = 'wide')
#EC = W_matrix(env.data = env.data,var.id = 'T2M')
EC = W_matrix(env.data = env.data,var.id = var)
distances = env_kernel(env.data = ET,gaussian = T)[[2]] 
kinship   = env_kernel(env.data = EC,gaussian = F, sd.tol = 3)[[2]]

# plot

require(superheat)
require(viridis)
## other plot
x11()
superheat(distances,
          pretty.order.rows = TRUE,
          pretty.order.cols = FALSE,
          row.dendrogram = T,
          col.dendrogram = T,
          grid.vline.col = "white",
          grid.hline.col = "white",
          #row.dendrogram = T,
          legend.width = 4,
          left.label.size = 0.1,
          bottom.label.text.size = 5,
          bottom.label.size = 0.2,
          bottom.label.text.angle = 90,
          heat.pal = viridis::inferno(100),
          #heat.pal = viridis::magma(100),
          legend.text.size = 17,
          #   X.text = round(as.matrix(a),1),X.text.col="white",
          legend.height=0.2)

var = colnames(env.data[c(9:17)])
