####### Packages #######
packages<- c(
  "Matrix", "fda", "fda.usc", "gstat", "splines", "colorspace", "sp", "geoR", "ggplot2",
  "reshape2", "readr", "readxl", "ggpubr", "ggrepel", "gridExtra", "patchwork", "purrr",
  "tidyr", "dplyr", "ggpubr", "rainbow", "cowplot", "readxl", "tidyverse", "grid", "lattice",
  "roahd", "ftsa", "fields", "SpatFD", "gstat", "sp", "sf", "reshape", "plotly", "processx"
)


for(p in packages){
  if(!require(p,character.only = TRUE)) install.packages(p, dependencies = TRUE)
  sapply(p, require, character.only = TRUE)
}

####### Working Directory #######
setwd("D:/WWU/Thesis")
source("Package.R")
setwd("D:/WWU/Thesis/Data/Dataset")

### Read dataset
for (i in 1:23) {
  if(i!=10){
    assign(paste0("s",i),value = read_excel(paste0("sujetoEpsd",i,".xlsx")))
    assign(paste0("s",i),value=as.matrix(dplyr::select(get(paste0("s",i)),-FRECUENCIA),row.names=1))
  }
}

### Create parameters and names for the data.
p = 228 ; nelec = 21 ; nvow = 5
vowels = c("a","e","i","o","u")

### Creation of objects for data functions. Given the functions created for the analysis for each vowel.
for (i in 1:23) {
  if(i!=10){
    assign(paste0("s",i,".list"),value= data.list(get(paste0("s",i)),p))
    assign(paste0("s",i,".gfdata"),value= gfdata(get(paste0("s",i)),p))
    assign(paste0("names(s",i,".list)"),value= vowels)
    assign(paste0("names(s",i,".gfdata)"),value= get(paste0("names(s",i,".list)")))
  }
}

### RMSE and MSE
nb = seq(10,30,1)
for (i in 1:23) {
  if(i!=10){
    assign(paste0("nb.MSE",i),value= gf.RMSE(get(paste0("s",i,".gfdata")),get(paste0("s",i,".list")),seq.nb = nb))
    #assign(paste0("nb.RMSE",i),value= gf.RMSE(get(paste0("s",i,".gfdata")),get(paste0("s",i,".list")),seq.nb = nb,option = "RMSE"))
  }
}

### Plot RMSE n basis selection
x11()
ggplot.RMSE(nb.MSE4,nb)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+ theme_bw()

n.basis<-c(11,13,12,13,11)
for (i in 1:23) {
  if(i!=10){
    assign(paste0("s",i,".gfd"),value=gfd.individual(get(paste0("s",i,".gfdata")),n.basis))}
}

### Change names of the data.
for (i in 1:23) {
  if(i!=10){
    nam<-paste0("names(s",i,".gfd)<-vowels")
    eval(parse(text = nam))
  }
}

plot<-NULL
for (i in 1:length(vowels)) {
    assign(paste0("plots",i),value=ggplot.gfd(s4.gfd,vow = vowels[i],trial = 3,legend = T))
  plot<-paste0(plot,"plots",i,",")
}

### Plots
plots<-paste0("plots<-ggarrange(",plot,"ncol=",3,",nrow=",2,")")
eval(parse(text = plots))
{x11()
  plots}

####### Descriptive statistics #######

### Calculate summary statistics
summary.s4.gfd<-summary.gfd(list.gfd=s4.gfdata)

### Plots mean and variance
plot<-NULL
for (i in 1:length(vowels)) {
  assign(paste0("plots",i),value=ggplot_combinedmean(s4.gfd,summary.s4.gfd,vow = vowels[i],trial = 3))
  plot<-paste0(plot,"plots",i,",")
}

plots<-paste0("plots<-ggarrange(",plot,"ncol=",3,",nrow=",2,")")
eval(parse(text = plots))
{x11()
  plots}

### Plots trimmed mean and median
plot<-NULL
for (i in 1:length(vowels)) {
  assign(paste0("plots",i),value=ggplot.gfdata.trim(summary.s4.gfd,vow = vowels[[i]],trial = 3, legend= T))
  plot<-paste0(plot,"plots",i,",")
}

plots<-paste0("plots<-ggarrange(",plot,"ncol=",3,",nrow=",2,")")
eval(parse(text = plots))
{x11()
  plots}

### Correlation coefficient 
for (i in 1:length(vowels)) {
  assign(paste0("s4.cor.",vowels[i]),value=cor.fd(1:p,get(paste0("s4.gfd"))[[vowels[i]]][['3']]))
}

windows()
par(mfrow=c(2,3))
for (i in 1:length(vowels)) {
  image.plot(1:p, 1:p, get(paste0("s4.cor.",vowels[i])), xlab='Hz', ylab='Hz')
}

### Outliers
for (i in 1:length(vowels)) {
  assign(paste0("fds.s4.",vowels[i]),value=fds(1:p,get(paste0("s4.list"))[[vowels[i]]][['3']]))
}

### Plot functional boxplot type HDR
windows()
par(mfrow=c(2,3))
for (i in 1:length(vowels)) {
  fboxplot( get(paste0("fds.s4.",vowels[i])), plot.type = "functional", type ="hdr", projmethod = "PCAproj",ylab="Brain Signal",xlab="Hz",na.rm = TRUE)
}

### Plot functional boxplot type bag
windows()
par(mfrow=c(2,3))
for (i in 1:length(vowels)) {
  fboxplot( get(paste0("fds.s4.",vowels[i])), plot.type = "functional", type ="bag", projmethod = "PCAproj",ylab="Brain Signal",xlab="Hz",na.rm = TRUE)
}


####### FPCA #######

s4.fpca =gfd_pca(s4.gfd)

### Plots PC functions
windows()
par(mfrow=c(2,3))
for (i in 1:length(vowels)) {
  plot(get(paste0("s4.fpca"))[[vowels[i]]][['3']][["harmonics"]], lwd = 1,xlab= "Hz", ylab='PCs')
}

### Plot scores
electrodes=c("1","2","3","4","5","6","7","8","9",
             "10","11","12","13","14","15","16","17",
             "18","19","20","21")

plot<-NULL
for (i in 1:length(vowels)) {
  assign(paste0("plots",i),value=ggplot.fpca.scores(s4.fpca,vow = vowels[i],trial=3,electrodes))
  plot<-paste0(plot,"plots",i,",")
}

plots<-paste0("plots<-ggarrange(",plot,"ncol=",3,",nrow=",2,")")
eval(parse(text = plots))
{x11()
  plots}

####### Functional geostatistics #######

coordenadas <- as.data.frame(read_table2("coordenadas.txt"))

for (i in 1:length(vowels)) {
  assign(paste0("vow",vowels[i],"_rep3"),value=t(get(paste0("s4.gfdata"))[[vowels[i]]][['3']]$data))
  assign(paste0("vow",vowels[i],"_rep3"),value=as.data.frame(get(paste0("vow",vowels[i],"_rep"))))
}

### Spatial functional data construction
for (i in 1:length(vowels)) {
  assign(paste0("sfd_",vowels[i]),value=SpatFD(get(paste0("vow",vowels[i],"_rep3"))%>% select(`E1`:`E21`),
                                              coords = coordenadas, basis = "Bsplines", nbasis = n.basis[i],
                                              lambda = 0.00001, nharm = 2,name =  paste0(i)))
  assign(paste0("ptj_",vowels[i]),value=data.frame(scores(get(paste0("sfd_",vowels[i])))[[1]]))
}

colnames(ptj_a)=colnames(ptj_e)=colnames(ptj_i)=colnames(ptj_o)=colnames(ptj_u)=c("x","y","f1","f2")
coordinates(ptj_a)=coordinates(ptj_e)=coordinates(ptj_i)=coordinates(ptj_o)=coordinates(ptj_u)=c("x","y")

### Variograms
var_af11=variogram(f1~1,ptj_a,cutoff=50,cressie=T)
var_af21=variogram(f2~1,ptj_a,cutoff=60,cressie=T)

var_ef11=variogram(f1~1,ptj_e,cutoff=50,cressie=T)
var_ef21=variogram(f2~1,ptj_e,cutoff=50,cressie=T)

var_if11=variogram(f1~1,ptj_i,cutoff=50,cressie=T)
var_if21=variogram(f2~1,ptj_i,cutoff=50,cressie=T)

var_of11=variogram(f1~1,ptj_o,cutoff=60,cressie=T)
var_of21=variogram(f2~1,ptj_o,cutoff=60,cressie=T)

var_uf11=variogram(f1~1,ptj_u,cutoff=60,cressie=T)
var_uf21=variogram(f2~1,ptj_u,cutoff=50,cressie=T)

### Models
models_a<- list(vgm(105, "Hol", 12),vgm(44.85, "Hol", 16.08))
models_e<- list(vgm(59.89, "Wav", 34.08),vgm(11.99, "Wav", 49.08))
models_i<- list(vgm(55.6, "Wav", 32.24),vgm(3.18, "Wav", 31.17))
models_o<- list(vgm(10.53, "Per", 154.67),vgm(4.98, "Hol", 11.60))
models_u<- list(vgm(19.92, "Per",161.17),vgm(4.55, "Per", 161.17))

### Spatial functional prediction
for (i in 1:length(vowels)) {
  assign(paste0("KS_SFD_",vowels[i],"_l"),value=KS_scores_lambdas(get(paste0("sfd_",vowels[i])), coordenadas ,
                                                                  method = "lambda", model = get(paste0("models_",vowels[i]))))
}

class(KS_SFD_a_l)
summary(KS_SFD_a_l)
recons_fd(KS_SFD_a_l)

for (i in 1:length(vowels)) {
  setwd(paste0("D:/WWU/Thesis/pred/",vowels[i]))
  crossval_loo(paste0("KS_SFD_",vowels[i],"_l"))
}


# Cross Validation
setwd("D:/WWU/Thesis/pred/a")
crossval_loo(KS_SFD_a_l)

setwd("D:/WWU/Thesis/pred/e")
crossval_loo(KS_SFD_e_l)

setwd("D:/WWU/Thesis/pred/i")
crossval_loo(KS_SFD_i_l)

setwd("D:/WWU/Thesis/pred/o")
crossval_loo(KS_SFD_o_l)

setwd("D:/WWU/Thesis/pred/u")
crossval_loo(KS_SFD_u_l)

### Kriging maps
setwd("D:/WWU/Thesis/krig")

dir.create(dirname(paste0('./a/')), recursive=TRUE)
dir.create(dirname(paste0('./e/')), recursive=TRUE)
dir.create(dirname(paste0('./i/')), recursive=TRUE)
dir.create(dirname(paste0('./o/')), recursive=TRUE)
dir.create(dirname(paste0('./u/')), recursive=TRUE)


Sys.setenv(PATH=paste0("C:/Users/bibia/AppData/Local/Programs/orca", Sys.getenv("PATH")))

N_img = 228

t = seq(sfd_a[[1]]$data_fd$basis$rangeval[1],
        sfd_a[[1]]$data_fd$basis$rangeval[2],length.out = N_img)


for (i in 1:N_img){
  orca(ggmap_KS(KS_SFD_a_l,
                map_path = NULL,
                window_time = t[i],
                method = "lambda")[[1]],  paste0('./a/vowa_',i,'.png'))
}

t = seq(sfd_e[[1]]$data_fd$basis$rangeval[1],
        sfd_e[[1]]$data_fd$basis$rangeval[2],length.out = N_img)

for (i in 1:N_img){
  orca(ggmap_KS(KS_SFD_e_l,
                map_path = NULL,
                window_time = t[i],
                method = "lambda")[[1]],  paste0('./e/vowe_',i,'.png'))
}

t = seq(sfd_i[[1]]$data_fd$basis$rangeval[1],
        sfd_i[[1]]$data_fd$basis$rangeval[2],length.out = N_img)

for (i in 1:N_img){
  orca(ggmap_KS(KS_SFD_i_l,
                map_path = NULL,
                window_time = t[i],
                method = "lambda")[[1]],  paste0('./i/vowi_',i,'.png'))
}

t = seq(sfd_o[[1]]$data_fd$basis$rangeval[1],
        sfd_o[[1]]$data_fd$basis$rangeval[2],length.out = N_img)

for (i in 1:N_img){
  orca(ggmap_KS(KS_SFD_o_l,
                map_path = NULL,
                window_time = t[i],
                method = "lambda")[[1]],  paste0('./o/vowo_',i,'.png'))
}


t = seq(sfd_u[[1]]$data_fd$basis$rangeval[1],
        sfd_u[[1]]$data_fd$basis$rangeval[2],length.out = N_img)

for (i in 1:N_img){
  orca(ggmap_KS(KS_SFD_u_l,
                map_path = NULL,
                window_time = t[i],
                method = "lambda")[[1]],  paste0('./u/vowu_',i,'.png'))
}

