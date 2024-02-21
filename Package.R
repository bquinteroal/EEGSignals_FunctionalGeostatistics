library(pracma)
#############################################
#######       Additional fuctions     #######
#############################################
#
##### Evaluates whether a matrix is positive definite
### matrix is a square matrix
is_positive_definite<-function(matrix){
  eigen<-eigen(matrix)$values
  con<-0
  for (i in 1:length(eigen)) {
    if(eigen[i]>0){con<-con+1}else{con<-con+0}
  }
  if(con==length(eigen)){return(TRUE)}else{FALSE}
}
#############################################
#######        Functional data        #######
#############################################
#
##### Create a list of the vowels and each of their repetitions.
#### "data": is a matrix object.
#### "p": is a scalar for dividing the data.
### The last column of the data matrix is the ranking according to the vowel of each row.
data.list=function(data,p){
  classes = levels(as.factor(data[,ncol(data)]))
  data.orig = list()
  for (i in 1:length(classes)) {
    data.orig[[classes[i]]]=subset(data[,-ncol(data)],data[,ncol(data)]==classes[i])
    n=NULL
    n[i]=nrow(data.orig[[i]])/p
    data.orig[[classes[i]]]=split.data.frame(data.orig[[i]], rep(1:n[i],each=p))
  }
  return( data.orig)
}
#
##### Construction of the data in a gfdata object of the fda.usc package according to each vowel and each repetition.
##### It is necessary to build the functional data with the fda.usc package.
#### "data": is a matrix object.
#### "p": is a scalar for dividing the data.
gfdata=function(data,p){
  if(!inherits(data, "matrix"))
    stop("'data' is not of class 'matrix'")
  #call=match.call()
  classes=levels(as.factor(data[,ncol(data)]))
  class=list()
  for (i in 1:length(classes)) {
    class[[i]]=subset(data[,-ncol(data)],data[,ncol(data)]==classes[i])
    if((nrow(class[[i]])/p)%%1!=0)
      stop(paste(" 'p' is not multiple of the number of rows in class",i))
    n=NULL
    n[i]=nrow(class[[i]])/p
    class[[i]]=split.data.frame(class[[i]], rep(1:n[i],each=p))
    class[[i]]=lapply(class[[i]], t)
    class[[i]]=lapply(class[[i]], fdata)
  }
  names(class)=paste("Class ", unique(as.factor(data[,ncol(data)])))
  class(class)="gfdata"
  return(class)
}
#
##### Construction of functional data
#### "gfdata": Created in function gfdata it must be this type of object.
#### "n.basis": Number of basis which is used in create.basis function of package fda.usc
#### "basis.type": Type of basis. A function create."type.basis".basis must exists. By default, bspline basis is used.
gfd=function(gfdata, n.basis, basis.type="bspline"){
  # if(!inherits(gfdata, "gfdata"))
  #   stop("'data' is not of class 'gfdata'")
  gfd=list()
  for(i in 1:length(gfdata)){
    gfd[[i]]=lapply(gfdata[[i]], fdata2fd, nbasis=n.basis, type.basis=basis.type)
  }
  class(gfd)="gfd"
  return(gfd)
}
#
#
##### Construction of functional data
#### "gfdata": Created in function gfdata it must be this type of object.
#### "n.basis": Number of basis which is used in create.basis function of package fda.usc
#### "basis.type": Type of basis. A function create."type.basis".basis must exists. By default, bspline basis is used.
gfd.individual=function(gfdata, n.basis, basis.type="bspline"){
  # if(!inherits(gfdata, "gfdata"))
  #   stop("'data' is not of class 'gfdata'")
  gfd=list()
  for(i in 1:length(gfdata)){
    gfd[[i]]=lapply(gfdata[[i]], fdata2fd, nbasis=n.basis[i], type.basis=basis.type)
  }
  class(gfd)="gfd"
  return(gfd)
}
#
##### Generates a summary list of the mean and variance of functional data
#### "list.gfd": object generated in the function "gfd".
### List Summary object gfd
### Generates a list with mean and variance for each functional data set.
summary.gfd=function(list.gfd){
  summary_gfd<-list()
  mean.gfd<-list()
  for (i in 1:length(list.gfd)) {
    mean.gfd[[i]]<-lapply(list.gfd[[i]], func.mean)
  }
  names(mean.gfd)<-names(list.gfd)
  
  var.gfd<-list()
  for (i in 1:length(list.gfd)) {
    var.gfd[[i]]<-lapply(list.gfd[[i]], func.var)
  }
  names(var.gfd)<-names(list.gfd)
  
  trim.mode.gfd<-list()
  for (i in 1:length(list.gfd)) {
    trim.mode.gfd[[i]]<-lapply(list.gfd[[i]], func.trim.mode)
  }
  names(trim.mode.gfd)<-names(list.gfd)
  
  
  trim.RP.gfd<-list()
  for (i in 1:length(list.gfd)) {
    trim.RP.gfd[[i]]<-lapply(list.gfd[[i]], func.trim.RP)
  }
  names(trim.RP.gfd)<-names(list.gfd)
  
  med.mode.gfd<-list()
  for (i in 1:length(list.gfd)) {
    med.mode.gfd[[i]]<-lapply(list.gfd[[i]], func.med.mode)
  }
  names(med.mode.gfd)<-names(list.gfd)
  
  
  med.RP.gfd<-list()
  for (i in 1:length(list.gfd)) {
    med.RP.gfd[[i]]<-lapply(list.gfd[[i]], func.med.RP)
  }
  names(med.RP.gfd)<-names(list.gfd)
  
  summary_gfd[[1]]<-mean.gfd
  summary_gfd[[2]]<-var.gfd
  summary_gfd[[3]]<-trim.mode.gfd
  summary_gfd[[4]]<-trim.RP.gfd
  summary_gfd[[5]]<-med.mode.gfd
  summary_gfd[[6]]<-med.RP.gfd
  names(summary_gfd)<-c("Mean","Var","Trim.mode","Trim.RP","Med.mode","Med.RP")
  return(summary_gfd)
}
#
##### Ouliers object gfd
#### "data.gfd": object generated in the function "gfd".
#### "vowel": vowel of interest in the data.gfd
#### Generates a list with the functional atypicals of each trial
outliers.gfd<-function(data.gfd,vowel){
  data.gfd.vowel<-data.gfd[[vowel]]
  outliers<-lapply(data.gfd.vowel,function(x){Out<-outliers.depth.trim(x)
  return(Out$outliers)})
  return(outliers)
}
#
##### Evaluate a list functional data object at specified argument values.
#### "grid": a vector or matrix of argument values at which the functional data object is to be evaluated.
#### "gfd_data": object generated in the function "gfd". Object functional data.
gfd_eval=function(grid,gfd_data){
  #if(!inherits(data, "gfdata"))
  #stop("'data' is not of class 'gfdata'")
  gfd_eval=list()
  for(i in 1:length(gfd_data)){
    gfd_eval[[i]]=lapply(X = gfd_data[[i]], FUN = eval.fd,evalarg = grid)
  }
  return(gfd_eval)
}
#
##### Sampling a gfd object
#### "gfd_data": object generated in the function "gfd". Object functional data.
#### "prop.train": Proportion to be taken from the sample.
gfd_sample=function(gfd_data, prop.train){
  gfd_sample=list()
  for(i in 1:length(gfd_data)){
    ntrial=length(gfd_data[[i]])
    sample=sort(sample(seq(1,ntrial),round(ntrial*prop.train)))
    gfd_sample[[i]]=gfd_data[[i]][sample]
  }
  return(gfd_sample)
}
##### Object with training and test data
#### "gfd_data": object generated in the function "gfd". Object functional data.
#### "prop.train": Proportion to be taken from the sample.
gfd_clasif_data=function(gfd_data, prop.train){
  data=list()
  gfd_train=list()
  gfd_test=list()
  for(i in 1:length(gfd_data)){
    ntrial=length(gfd_data[[i]])
    sample=sort(sample(seq(1,ntrial),round(ntrial*prop.train)))
    gfd_train[[i]]=gfd_data[[i]][sample]
    gfd_test[[i]]=gfd_data[[i]][-sample]
  }
  data=list(train=gfd_train, test=gfd_test)
  return(data)
}
#
##### Sum of squares for the creation of functional data.
#### "gfdata": Created in function gfdata it must be this type of object.
#### "data": Original data
#### "nbasis": Number of basis which is used in create.basis function of package fda.usc
#### "btype": Type of basis. A function create."type.basis".basis must exists. By default, bspline basis is used.
gf_sse=function(gfdata, data, btype = "bspline",nbasis){
  gfd=gfd(gfdata,basis.type = btype, n.basis =nbasis)
  arvals=gfdata[[1]][[1]]$argvals
  eval=gfd_eval(arvals,gfd)
  MSE=matrix(NA,ncol = length(eval),nrow = length(eval[[1]]))
  RMSE=matrix(NA,ncol = length(eval),nrow = length(eval[[1]]))
  for (j in 1:length(eval)) {
    for(i in 1:length(eval[[1]])){
      MSE[i,j]=mean(apply((eval[[j]][[i]]-data[[j]][[i]])^2,2,mean))
      RMSE[i,j]=sqrt(MSE[i,j])
    }
  }
  return(list(MSE=MSE,RMSE=RMSE))
}
##### Sum of squares of the error for the creation of functional data.
#### "gfdata": Created in function gfdata it must be this type of object.
#### "data": Original data.
#### "seq.nb": Vector with the number of bases.
#### "option": Default "MSE".

gf.RMSE=function(gfdata,data,seq.nb, option = "MSE"){
  nb.sse = list()
  opt = list()
  for (i in 1:length(seq.nb)) {
    nb.sse[[i]] = gf_sse(gfdata = gfdata, data = data,
                         btype = "bspline",nbasis = seq.nb[i])
    opt[[i]] = nb.sse[[i]][[option]]
  }
  opt.mean = matrix(unlist(lapply(opt,apply,2,mean)),ncol = length(gfdata ),
                    nrow = length(seq.nb),byrow = T)
  colnames(opt.mean) = names(gfdata)
  return(opt.mean)
}
#
##### Calculates the average of the averages of each performance per vowel.
#### "data.train.pca": is a list of lists. Where each sub-list contains the FACP of each realization of each group..
mean.mean=function(data.train.pca){
  mean.vowels<-list()
  for (i in 1:length(data.train.pca)) {
    data.mean<-c()
    for (j in 1:length(data.train.pca[[i]])) {
      mean.vow<-data.train.pca[[i]][[j]]$meanfd
      eval<-eval.fd(seq(1,mean.vow$basis$rangeval[2]),mean.vow)
      data.mean<-rbind(data.mean,t(eval))
    }
    mean.gfd<-fdata(data.mean)
    mean.v<-func.mean(mean.gfd)
    mean.vowels[[i]]<-mean.v
  }
  names(mean.vowels)<-names(data.train.pca)
  return(mean.vowels)
}
#
#############################################
#######            SPATIAL            #######
#############################################
#
##### Given a sample of functional data generates functional principal components analysis.
#### "gfd_sample": is a list object of functional data.
### The last column of the data matrix is the ranking according to the vowel of each row.
gfd_pca=function(gfd_data){
  gfd_pca=map(gfd_data,function(vowel){map(vowel,pca.fd,centerfns=T)})
  return(gfd_pca)
}

plot.gfd_pca<-function(data,p,vowel,nco,nrow){
  data.list<-data.list(data=data,p=p)
  x11()
  par(mfrow=c(nrow,nco),mar=c(2,2,2,2))
  for (i in 1:length(data.list[[vowel]])) {
    fds.data<-fds(x=1:p,y=data.list[[vowel]][[i]])
    fboxplot(data=fds.data, plot.type = "functional", type = type, projmethod = "PCAproj",ylab="",xlab="",na.rm = TRUE)
  }
}

#
##### Functions to take objects generated by gfd_pca
#### harmonics a functional data object for the harmonics or eigenfunctions
#### values the complete set of eigenvalues
#### scores s matrix of scores on the principal components or harmonics
#### varprop a vector giving the proportion of variance explained by each eigenfunction
#### meanfd a functional data object giving the mean function
take_scores=function(gfd_pca) {ret=gfd_pca$scores; return(ret)}
#
take_scores1=function(gfd_pca) {ret=gfd_pca$scores[,1]; return(ret)}
take_scores2=function(gfd_pca) {ret=gfd_pca$scores[,2]; return(ret)}
#
take_varprop=function(gfd_pca) {ret=gfd_pca$varprop; return(ret)}
#
take_values=function(gfd_pca){ret1=gfd_pca$values;return(ret1)}
#
take_values1=function(gfd_pca){ret1=gfd_pca$values[[1]];return(ret1)}
#
take_values2=function(gfd_pca){ret2=gfd_pca$values[[2]];return(ret2)}
#
take_coords=function(ret, coord){
  data=data.frame(coord,ret)
  coordinates(data)=c("x","y")
  return(data)}

#
##### Create variogram of the geoR package for each repetition of each vowel and group in a list.
#### "gfd_pca_Data": is a list of lists. Where each sub-list contains the FACP of each realization of each group.
#### "coord": an n \times 2n×2 matrix containing coordinates of the nn data locations in each row. Defaults to geodata$coords, if provided.
#### "pairsmin": a integer number defining the minimum numbers of pairs for the bins. For option = "bin",
####             bins with number of pairs smaller than this value are ignored. Defaults to 2.
gfd_variog_geoR=function(gfd_pca_Data, coord, pairsmin=2){
  scores=list()
  dataframe=list()
  geodata.f1=list()
  geodata.f2=list()
  variog.f1=list()
  variog.f2=list()
  variogs=list()
  for(i in 1:length(gfd_pca_Data)){
    scores[[i]]=lapply(X=gfd_pca_Data[[i]], FUN=take_scores)
    dataframe[[i]]=lapply(scores[[i]],FUN=data.frame, coord)
    geodata.f1[[i]]=lapply(X=dataframe[[i]],FUN=as.geodata,coords.col=3:4, data.col=1)
    geodata.f2[[i]]=lapply(X=dataframe[[i]],FUN=as.geodata,coords.col=3:4, data.col=2)
    variog.f1[[i]]=lapply(X=geodata.f1[[i]], FUN=variog, pairs.min=pairsmin)
    variog.f2[[i]]=lapply(X=geodata.f2[[i]], FUN=variog, pairs.min=pairsmin)
    variogs[[i]]=list(variog.f1[[i]], variog.f2[[i]])
  }
  return(list(geodata1= geodata.f1,variogr=variog.f1, geodata2= geodata.f2,variogr2=variog.f2))
}
#
##### Create variogram of the gstat package for each repetition of each vowel and group in a list.
#### "gfd_pca_Data": is a list of lists. Where each sub-list contains the FACP of each realization of each group.
#### "coord": spatial data locations..
#### "model": model type, e.g. "Exp", "Sph", "Gau", "Mat". Calling vgm() without a model argument returns a data.frame with available models.
gfd_variog=function(gfd_pca_Data, coords, model){
  distance=as.matrix(dist(coords))
  scores=list()
  datageo=list()
  geodata.f1=list()
  geodata.f2=list()
  variog.f1=list()
  variog.f2=list()
  values.f1=list()
  values.f2=list()
  vgm.f1=list()
  vgm.f2=list()
  for(i in 1:length(gfd_pca_Data)){
    scores[[i]]=lapply(X=gfd_pca_Data[[i]], FUN=take_scores)
    datageo[[i]]=lapply(X=scores[[i]], FUN=take_coords, coord=coords)
    geodata.f1[[i]]=lapply(X=datageo[[i]],FUN=gstatF,id="X1", formula=X1~1)
    geodata.f2[[i]]=lapply(X=datageo[[i]],FUN=gstatF,id="X2", formula=X2~1)
    variog.f1[[i]]=lapply(X=geodata.f1[[i]], FUN=variogram, cutoff=max(distance))
    variog.f2[[i]]=lapply(X=geodata.f2[[i]], FUN=variogram, cutoff=max(distance))
    values.f1[[i]]=lapply(X=gfd_pca_Data[[i]], FUN=take_values1)
    values.f2[[i]]=lapply(X=gfd_pca_Data[[i]], FUN=take_values2)
    vgm.f1[[i]]=lapply(X=values.f1[[i]],FUN=vgm,model=model,range=max(distance)/4)
    vgm.f2[[i]]=lapply(X=values.f2[[i]],FUN=vgm,model=model,range=max(distance)/4)
    spatial=list(variog.f1=variog.f1, vgm.f1=vgm.f1,data=datageo)
  }
  return(spatial)
}
#
##### Create list for each vowel, Fit ranges and/or sills from a simple or nested variogram model to a sample variogram.
#### "variog": sample variogram, output of variogram, of the function "gfd_variog".
#### "vgm": variogram model, output of vgm; see Details below for details on how NA values in model are initialised.
gfd_fit=function(variog,vgm){
  fit.f1=list()
  for(i in 1:length(vgm)){
      fit.f1[[i]]=fit.variogram(variog[[i]], model=vgm[[i]],fit.sills =T)
  }
  return(fit.f1)
}
#
##### Create list for each vowel, Generates a semivariance values given a variogram model.
#### "gfd_variog": sample variogram, output of variogram, of the function "gfd_variog".
#### "coords": spatial data locations.
gfd_variogramLine=function(gfd_variog, coords){
  lines=list()
  for(i in 1:length(gfd_variog)){
    lines[[i]]=lapply(gfd_variog[[i]], FUN=variogramLine,dist_vector=as.matrix(dist(coords)), covariance=T)
  }
  return(lines)
}


#########################################
#######           plots           #######
#########################################
##### ggplot mean and var
#### "gfdata": list with mean and variance for each functional data set,function summary.gfd
#### "vow": Value of the vowel
#### "trial": Value of trial
#### "legend": if TRUE generates legend for each curve. default FALSE

ggplot.gfdata.mean = function(gfdata,vow,trial,legend=FALSE){
  eval=rbind(gfdata[["Mean"]][[vow]][[trial]]$data,gfdata[["Var"]][[vow]][[trial]]$data)
  rownames(eval)<-c("Mean","Var")
  melt_s=melt(eval);class(melt_s)= "data.frame"
  as.factor(melt_s$Var2)
  names(melt_s)=c("Summary","t","X_t")
  ggplot(melt_s,aes(x=t,y=X_t,col=Summary)) +
    geom_line()+labs(title=paste("Vowel ",vow))+
    xlab("Hz") +
    ylab("Brain Signal")+
    if(legend==FALSE){theme(legend.position = 'none')}
}


ggplot_combinedmean <- function(gfdata, sum,vow, trial, legend = FALSE) {
  # Original data
  arvals <- gfdata[[vow]][[trial]]$basis$rangeval[1]:gfdata[[vow]][[trial]]$basis$rangeval[2]
  eval_original <- eval.fd(arvals, gfdata[[vow]][[trial]])
  melt_original <- melt(eval_original)
  class(melt_original) <- "data.frame"
  as.factor(melt_original$Var2)
  names(melt_original) <- c("t", "Electrode", "X_t")
  
  # Mean and variance data
  eval_mean <- rbind(sum[["Mean"]][[vow]][[trial]]$data, sum[["Var"]][[vow]][[trial]]$data)
  rownames(eval_mean) <- c("Mean", "Var")
  melt_mean <- melt(eval_mean, as.is = TRUE)
  class(melt_mean) <- "data.frame"
  as.factor(melt_mean$Var2)
  names(melt_mean) <- c("Summary", "t", "X_t")
  
  # Plotting
  p <- ggplot() +
    geom_line(data = melt_original, aes(x = t, y = X_t, col = Electrode, group = Electrode), linewidth = 0.7) +
    geom_line(data = melt_mean, aes(x = t, y = X_t, col = Summary),linewidth = 1) +
    labs(title = paste("Vowel ", vow), x = "Hz", y = "Brain Signal") +
    scale_color_manual(values = c("gray","gray","gray","gray","gray","gray","gray","gray","gray",
                                  "gray","gray","gray","gray","gray","gray","gray","gray","gray",
                                  "gray","gray","gray","black", "red"), guide = FALSE) +
    theme_minimal() +
    theme(
      legend.position = ifelse(legend, "right", "none"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  return(p)
}


ggplot.gfdata.trim = function(gfdata,vow,trial,legend=FALSE){
  eval=rbind(gfdata[["Trim.mode"]][[vow]][[trial]]$data, gfdata[["Trim.RP"]][[vow]][[trial]]$data, gfdata[["Med.mode"]][[vow]][[trial]]$data,gfdata[["Med.RP"]][[vow]][[trial]]$data)
  rownames(eval)<-c("Trimmed mean (mode)","Trimmed mean (RP)","Median (mode)","Median (RP)")
  melt_s=melt(eval);class(melt_s)= "data.frame"
  as.factor(melt_s$Var4)
  names(melt_s)=c("Summary","t","X_t")
  ggplot(melt_s,aes(x=t,y=X_t,col=Summary)) +
    geom_line(size=0.71)+labs(title=paste("Vowel ",vow))+
    xlab("Hz") +
    ylab("Brain Signal")+
    scale_color_brewer(palette = "RdYlBu")+
    theme_minimal() +
    theme(
      legend.position = "none",  # Hide the legend if legend = FALSE
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Black panel border
    )
}

##### Plot object gfd outlier
#### "data.gfd": object generated in the function "gfd".
#### "outlier.list": object generated in the function "outlier.gfd"
#### "vowel": vowel from which to generate outliers per trial
plot.outlier.gfd<-function(data.gfd,outlier.list,vowel){
  n<-0
  for(i in 1:length(outlier.list)){
    outlier<-outlier.list[[i]]
    if(!(length(outlier)==0)){
      n<-n+1
    }
  }
  plotdim<-ceiling(sqrt(n))
  windows()
  par(mfrow=c(plotdim,plotdim))

  for(i in 1:length(outlier.list)){
    outlier<-outlier.list[[i]]
    if(!(length(outlier)==0)){
      plot(data.gfd[[vowel]][[i]],col=1,main=paste("Outlier rep",i),xlab="",ylab="", xaxt="n")
      lines(data.gfd[[vowel]][[i]][outlier.list[[i]]],col=2,lwd=2,lty=1)
      legend("topright",outlier.list[[i]],col=2,lwd=2,cex = 1,lty=1)
    }
  }
}
#
##### plot functional boxplots with outliers for package fds
#### "data": is a matrix object.
#### "p": is a scalar for dividing the data.
#### "vowel": vowel from which to generate outliers per trial and boxplot
#### "nco": Number of columns for graph
#### "nrow": Number of rows for graph
#### "type": Type of boxplot. When type = "bag", a bagplot is provided. When type = "hdr", a HDR boxplot is provided.
plot.outlier.fds<-function(data,p,vowel,nco,nrow,type){
  data.list<-data.list(data=data,p=p)
  x11()
  par(mfrow=c(nrow,nco),mar=c(2.5,2.5,2.5,2.5))
  for (i in 1:length(data.list[[vowel]])) {
    fds.data<-fds(x=1:p,y=data.list[[vowel]][[i]])
    par(ask=TRUE)
    fboxplot(data=fds.data, plot.type = "functional", type = type, projmethod = "PCAproj",ylab="",xlab="",na.rm = TRUE)
  }
}
##### plot for object gfdata
#### "gfdata": object generate for function gfdata
#### "vow": Value of the vowel
#### "trial": Value of trial
ggplot.gfdata = function(gfdata,vow,trial){
  eval=t(gfdata[[vow]][[trial]]$data)
  melt_s=melt(eval);class(melt_s)= "data.frame"
  as.factor(melt_s$Var2)
  names(melt_s)=c("t","Electrode","X_t")
  ggplot(melt_s,aes(x=t,y=X_t,col=Electrode)) +
    geom_line()+labs(title=paste("Vowel ",vow, " Trial ",trial))+
                     xlab("Hz") +
                       ylab("Brain Signal")
}
#
##### plot for object gfd
#### "gfdata": object generate for function gfd
#### "vow": Value of the vowel
#### "trial": Value of trial
#### "legend": if TRUE generates legend for each curve. default FALSE
ggplot.gfd=function(gfdata,vow,trial,legend=FALSE){
  colores=c("darkblue","darkcyan","darkgoldenrod","darkgray","#581845","darkgreen","darkkhaki","darkmagenta","darkolivegreen","darkorange","darkorchid","darkred","darksalmon","darkseagreen","darkslateblue","darkslategray","#008080","darkturquoise",	
            "darkviolet","deeppink","deepskyblue")
  arvals=gfdata[[vow]][[trial]]$basis$rangeval[1]:gfdata[[vow]][[trial]]$basis$rangeval[2]
  eval=eval.fd(arvals, gfdata[[vow]][[trial]])
  melt_s=melt(eval);class(melt_s)= "data.frame"
  as.factor(melt_s$Var2)
  names(melt_s)=c("t","Electrode","X_t")
  ggplot(melt_s, aes(x = t, y = X_t, col = Electrode#, group = Electrode
                     )) +
    geom_line(linewidth=0.75) +
    #scale_x_continuous(breaks=seq(0, 15000, 2000))+
    labs(title = paste("Vowel ", vow)) +
    xlab("Hz") +
    ylab("Brain Signal") +
    scale_color_manual(values = colores) +
    #scale_color_brewer(palette = "Set2")+
    #scale_color_manual(values = rep("gray", length(unique(melt_s$Electrode))), guide = FALSE) +
    theme_minimal() +
    theme(
      legend.position = "right",  # Hide the legend if legend = FALSE
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Black panel border
    )
}
##### plot for object gf.RMSE
#### "seq.nb": Vector with the number of bases.
#### "ylab": name for each graph
ggplot.RMSE=function(gf.RMSE,seq.nb,ylab = "Mean MSE" ){
  base = list()
  g = list()
  for (i in 1:ncol(gf.RMSE)) {
    base[[i]]=data.frame(nb=seq.nb,MSE=gf.RMSE[,i])
    g[[i]] = ggplot(base[[i]], aes(x=nb, y=MSE)) +
      geom_line() +
      labs(title = " ", subtitle = paste("Vowel",colnames(gf.RMSE)[i]),
           x = "Basis Number",  y = ylab)+
      scale_x_discrete(limits=seq.nb) +
      geom_point(aes(col=MSE)) +
      theme_minimal() +
      theme(
        legend.position = "none",  # Hide the legend if legend = FALSE
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Black panel border
      )

  }
  n <- length(g)
  nrow <- floor(sqrt(n))
  return(do.call("grid.arrange", c(g, nrow=nrow)))
}
##### ggplot for varition explicate for scores
#### "mean.values": vector of values for each vowel
#### "nv": number of values
ggplot.eig = function(mean.values, nv){
  eig = lapply(mean.values, function(a) data.frame(nva = 1:nv, Prop = a[1:nv]))
  names(eig) = vowels
  p = list()
  for (i in 1:length(eig)) {
    p[[i]] = ggplot(eig[[i]],aes(x = nva, y = Prop)) + geom_line() +
      geom_point(aes(col=Prop)) +
      labs(title = " ", subtitle = paste("Vowel", names(eig)[i]),
           x = "Eigenvalues Number",
           y = "Mean explained variance (%) ",
           col = "%") +
      scale_x_continuous(breaks = seq(from = 0, to = nv, by = 1))+
      theme_minimal() +
      theme(
        legend.position = "right",  # Hide the legend if legend = FALSE
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Black panel border
      )

  }
  n <- length(p)
  nrow <- floor(sqrt(n))
  return(do.call("grid.arrange", c(p, nrow=nrow)))
}

##### ggplot for scores

ggplot.fpca.scores = function(fpcdata,vow,trial,electrodes,legend=FALSE){
  eval=rbind(fpcdata[[vow]][[trial]]$scores)
  eval=as.data.frame(eval)
  colnames(eval)<-c("PCS1","PCS2")
  class(eval)= "data.frame"
  rownames(eval)<-electrodes
  ggplot(eval,aes(x=PCS1,y=PCS2)) +
    geom_point(size=2)+geom_text(label=rownames(eval),nudge_x = 0.4, nudge_y = 0.4,size=3.2,
                                 check_overlap = T)+labs(title=paste("Vowel ",vow))+
    geom_hline(yintercept = 0, color = "gray")+ 
    geom_vline(aes(xintercept = 0), color = "gray")+ 
    xlab("PC Score 1") +
    ylab("PC Score 2")+
    theme_minimal() +
    theme(
      legend.position = "none",  # Hide the legend if legend = FALSE
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # Black panel border
    )
}


###################################################
#####       Distances between vectors     #########
###################################################
#
#### Calculates the Euclidean distance between vectors
### a is a vector of dimension n
### b is a vector of dimension n
euclidean_distance = function(a, b){
  # it is found that they have the same number of observations
  if(length(a) == length(b)){
    sqrt(sum((a-b)^2))
  } else{
    stop('Vectors must be of the same length')
  }
}
#
#
#### Calculates the mahattan distance between vectors
### a is a vector of dimension n
### b is a vector of dimension n
manhattan_distance = function(a, b){
  # It is verified that they have the same number of observations.
  if(length(a) == length(b)){
    sum(abs(a-b))
  } else{
    stop('Vectors must be of the same length')
  }
}
#
#
#### Calculates the similarity cosine distance between vectors
### a is a vector of dimension n
### b is a vector of dimension n
cos_similarity = function(a,b){
  if(length(a) == length(b)){
    num = sum(a *b, na.rm = T)
    den = sqrt(sum(a^2, na.rm = T)) * sqrt(sum(b^2, na.rm = T))
    result = num/den

    return(1-result)
  } else{
    stop('Vectors must be of the same length')
  }
}
#
#
#### Calculates the minkowski distance between vectors
### a is a vector of dimension n
### b is a vector of dimension n
### p  is a denominator in the exponent 1/p
minkowski_distance = function(a,b,p){
  if(p<=0){
    stop('p must be higher than 0')
  }

  if(length(a)== length(b)){
    sum(abs(a-b)^p)^(1/p)
  }else{
    stop('Vectors must be of the same length')

  }
}
#
#
#### Calculates the mahalanobis distance between vectors
### vector is a vector of dimension n
### compa is a vector of dimension n a compare
### cov matrix of covariance for the process of distance
mahalanobis.distance=function(vector,compa,cov){
  dist<-t(compa-vector)%*%Matrix::solve(cov)%*%(compa-vector)
  return(dist)
}

#########################################
#######      Classification       #######
#########################################
####

clasification<-function(data.train.pca,new.data,k,distance,mcov=NULL,n.basis, basis.type){

  mconf<-matrix(NA,nrow = length(data.train.pca[[1]]),ncol = length(data.train.pca))
  for (i in 1:length(data.train.pca)) {
    new.basis<-fdata2fd(fdata(t(new.data)), nbasis=n.basis[[i]], type.basis=basis.type)
    for (j in 1:length(data.train.pca[[1]])) {

      new.vector<-c()
      for (k in 1:dim(new.basis$coefs)[2]) {
        new.vector[k]<-inprod((new.basis[k]-data.train.pca[[i]][[j]]$meanfd),data.train.pca[[i]][[j]][[1]][1])
      }

      if(distance=="mahalanobis"){
        mconf[j,i]<-mahalanobis.distance(data.train.pca[[i]][[j]][["scores"]][,1],new.vector,mcov[[i]])
      }
      if(distance=="euclidean"){
        mconf[j,i]<-euclidean_distance(new.vector,data.train.pca[[i]][[j]][["scores"]][,1])
      }
      if(distance=="manhattan"){
        mconf[j,i]<-manhattan_distance(new.vector,data.train.pca[[i]][[j]][["scores"]][,1])
      }
      if(distance=="cos_similarity"){
        mconf[j,i]<-cos_similarity(new.vector,data.train.pca[[i]][[j]][["scores"]][,1])
      }
      if(distance=="minkowski"){
        mconf[j,i]<-minkowski_distance(new.vector,data.train.pca[[i]][[j]][["scores"]][,1],1)
      }

    }
  }

  mconf1<-as.data.frame(mconf)
  mconf2<-mconf1%>%pivot_longer(cols = starts_with("V"),names_to = "Vowel")
  mconf3<-mconf2[order(mconf2$value),]

  t<-table(mconf3[1:k,1])

  sort(as.vector(t/sum(t)),decreasing = TRUE)[1]
  clas0<-as.data.frame(t/sum(t))
  clas<-as.character(clas0[order(clas0$Freq,decreasing = TRUE),][1,1])
  return(clas)
}












