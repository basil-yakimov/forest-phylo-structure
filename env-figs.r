library(readxl)

env <- read_excel("raw.data/dls_envi_biomodel_transect.xlsx", sheet = 1)


plot(T ~ ele, data = env, pch = 21, bg = rgb(runif(1), runif(1), runif(1), runif(1)))
plot(WMT ~ ele, data = env, pch = 21, bg = rgb(runif(1), runif(1), runif(1), runif(1)))
plot(CMT ~ ele, data = env, pch = 21, bg = rgb(runif(1), runif(1), runif(1), runif(1)))
plot(P ~ ele, data = env, pch = 21, bg = rgb(runif(1), runif(1), runif(1), runif(1)))
plot(PE ~ ele, data = env, pch = 21, bg = rgb(runif(1), runif(1), runif(1), runif(1)))
