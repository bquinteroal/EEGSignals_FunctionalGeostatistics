crossval_loo = function(object,plot_show=TRUE){
  
  # validation -------------------------------------------------------------
  if(!inherits(object,"KS_pred")){
    stop("input object must be KS_pred")
  }
  
  # object ----------------------------------------------------------------
  SFD_main = object$SFD
  KS_main = object
  args_SFD = SFD_main[[1]]$call_args
  
  data = args_SFD$data
  data_names=names(data)
  coord = args_SFD$coords
  n = ncol(data) 
  
  residual_norm_lambda = rep(NA,n)
  residual_norm_scores = rep(NA,n)
  
  # iteration on method both ------------------------------------------------
  if (!is.null(KS_main$KS_scores)&&!is.null(KS_main$KS_lambda)){
    for (i in 1:n){
      
      new_coord_i = coord[i,]
      data_i = data[,-i]
      coord_i = coord[-i,]
      
      suppressWarnings({
      SFD_i = SpatFD(
        data = data_i,
        coords = coord_i,
        basis = args_SFD$basis,
        nbasis = args_SFD$nbasis,
        lambda = args_SFD$lambda,
        nharm = args_SFD$nharm,
        name = args_SFD$name,
        add = args_SFD$add)
      
      KS_i_lambda = KS_scores_lambdas(
        SFD_i, 
        new_coord_i, 
        method = "lambda", 
        model = object$model)
      KS_i_scores = KS_scores_lambdas(
        SFD_i, 
        new_coord_i, 
        method = "scores", 
        model = object$model)
      predict_i_lambda = recons_fd(KS_i_lambda)
      predict_i_scores = recons_fd(KS_i_scores)
      })
      
      residual_norm_lambda[i]=sqrt(inprod(SFD_main[[1]]$data_fd[i]-predict_i_lambda,SFD_main[[1]]$data_fd[i]-predict_i_lambda,rng=c(1,nrow(data))))
      residual_norm_scores[i]=sqrt(inprod(SFD_main[[1]]$data_fd[i]-predict_i_scores,SFD_main[[1]]$data_fd[i]-predict_i_scores,rng=c(1,nrow(data))))
      
      if(plot_show){
        par(ask=TRUE)
        plot(SFD_main[[1]]$data_fd[i], las=2)
        par(new=TRUE)
        plot(predict_i_lambda, ann=FALSE, axes=FALSE,col=2)
        par(new=TRUE)
        plot(predict_i_scores, ann=FALSE, axes=FALSE,col=3)
        par(new=FALSE)
        legend("topright", legend = c("KS lambda", "KS scores", data_names[i]), col = c("red","green" ,"black"), lty = 1)
        #par(ask=FALSE)
      }
    }
  }
  
  # iteration on method lambda ------------------------------------------------
  if (is.null(KS_main$KS_scores)&&!is.null(KS_main$KS_lambda)){
    
    for (i in 1:n){
      
      # Construct the filename based on the value of i
      pdf_filename <- paste("plot_", i, ".pdf")
      
      # Open a PDF file for the current plot
      pdf(pdf_filename)
      
      new_coord_i = coord[i,]
      data_i = data[,-i]
      coord_i = coord[-i,]
      
      suppressWarnings({
      SFD_i = SpatFD(
        data = data_i,
        coords = coord_i,
        basis = args_SFD$basis,
        nbasis = args_SFD$nbasis,
        lambda = args_SFD$lambda,
        nharm = args_SFD$nharm,
        name = args_SFD$name,
        add = args_SFD$add)
      
      KS_i_lambda = KS_scores_lambdas(
        SFD_i, 
        new_coord_i, 
        method = "lambda", 
        model = object$model)
      predict_i_lambda = recons_fd(KS_i_lambda)
      })
      
      residual_norm_lambda[i]=sqrt(inprod(SFD_main[[1]]$data_fd[i]-predict_i_lambda,SFD_main[[1]]$data_fd[i]-predict_i_lambda,rng=c(1,nrow(data))))
      
      if(plot_show){
        
        #par(ask=TRUE)
        par(mfrow = c(1,1))
        plot(SFD_main[[1]]$data_fd[i], las=2,xlab= " ", ylab=" ",  cex.lab=1.5, cex.axis=1.5)
        par(mfrow = c(1,1), new=TRUE)
        plot(predict_i_lambda, ann=FALSE, axes=FALSE,col=2)
        par(mfrow = c(1,1), new=FALSE)
        legend("topright", legend = c("Prediction", data_names[i]), col = c("red", "black"), lty = 1,pt.cex = 1)
        #par(ask=FALSE)
      }
      dev.off()
    }
  }
  
  # iteration on method scores ------------------------------------------------
  if (!is.null(KS_main$KS_scores)&&is.null(KS_main$KS_lambda)){
    for (i in 1:n){
      
      new_coord_i = coord[i,]
      data_i = data[,-i]
      coord_i = coord[-i,]
      
      suppressWarnings({
      SFD_i = SpatFD(
        data = data_i,
        coords = coord_i,
        basis = args_SFD$basis,
        nbasis = args_SFD$nbasis,
        lambda = args_SFD$lambda,
        nharm = args_SFD$nharm,
        name = args_SFD$name,
        add = args_SFD$add)
      
      KS_i_scores = KS_scores_lambdas(
        SFD_i, 
        new_coord_i, 
        method = "scores", 
        model = object$model)
      predict_i_scores = recons_fd(KS_i_scores)
      })
      
      residual_norm_scores[i]=sqrt(inprod(SFD_main[[1]]$data_fd[i]-predict_i_scores,SFD_main[[1]]$data_fd[i]-predict_i_scores,rng=c(1,nrow(data))))
      
      if(plot_show){
        par(ask=TRUE)
        par(mfrow = c(1,1))
        plot(SFD_main[[1]]$data_fd[i], las=2)
        par(mfrow = c(1,1), new=TRUE)
        plot(predict_i_scores, ann=FALSE, axes=FALSE,col=3)
        par(mfrow = c(1,1), new=FALSE)
        legend("topright", legend = c("KS scores", data_names[i]), col = c("green", "black"), lty = 1)
        par(ask=FALSE)
      }
    }
  }
  
  # output of mean residual norm ----------------------------------------------
  cat("## Mean Residual Norm","\n")
  cat("lambda method: ", mean(residual_norm_lambda),"<---\n")
  print(summary(residual_norm_lambda))
  print(residual_norm_lambda)
  cat("scores method: ", mean(residual_norm_scores),"<---\n")
  print(summary(residual_norm_scores))
}
