ggmap_KS <-
  function(KS, map_path=NULL, window_time = NULL, method = "lambda", map_n = 5000, zmin = NULL, zmax = NULL, graph = "plotly"){
    if (!is.null(map_path)){
      if(is.character(map_path)){
        map <- sf::st_read(map_path)
      }else{
        map = map_path
      }
    }else{
      mx <- min(KS$SFD[[1]]$coords[,1])
      Mx <- max(KS$SFD[[1]]$coords[,1])
      my <- min(KS$SFD[[1]]$coords[,2])
      My <- max(KS$SFD[[1]]$coords[,2])
      map <- sf::st_polygon(list(
        matrix(c(mx,my,Mx,my,Mx,My,mx,My,mx,my),byrow = T,ncol = 2)),
      )
    }
    
    newcoords <- sf::st_sample(map, map_n, type = "regular")
    newcoords <- sf::st_coordinates(newcoords)
    colnames(newcoords) <- colnames(KS$SFD[[1]]$coords)
    
    KS_SFD <- KS_scores_lambdas(KS$SFD, newcoords, model = KS$model, method = method, name = KS$name)
    
    SFD <- recons_fd(KS_SFD)
    
    if(is.null(window_time)) {
      times <- SFD$basis$rangeval[1]
    } else if (!(all(window_time >= SFD$basis$rangeval[1]) && all(window_time <= SFD$basis$rangeval[2]))) {
      stop(paste("window_time is out of bounds: Must be some value(s) between ", SFD$basis$rangeval[1], "and ",SFD$basis$rangeval[2]))
    }  else {
      times <- sort(window_time)
    }
    
    eval <- fda::eval.fd(times, SFD)
    
    melt_s <- suppressWarnings(reshape::melt(eval))
    
    melt_s$X2 <- as.factor(melt_s$X2)
    
    melt_s$X <- as.factor(melt_s$X2)
    levels(melt_s$X) <- newcoords[,1]
    
    melt_s$Y <- as.factor(melt_s$X2)
    levels(melt_s$Y) <- newcoords[,2]
    
    names(melt_s) = c("Time","Prediction","Value", "X", "Y")
    
    melt_s$Time <- as.factor(melt_s$Time)
    
    graf <- list()
    if(is.null(zmin)){zmin = min(melt_s$Value)}
    if(is.null(zmax)){zmax = max(melt_s$Value)}
    
    rownames(KS$SFD[[1]]$coords)<-c("E1","E2","E3","E4","E5","E6","E7","E8","E9",
                      "E10","E11","E12","E13","E14","E15","E16","E17",
                      "E18","E19","E20","E21")
    
    for(i in 1:nlevels(melt_s$Time)){
      
      melt_s_2 <- melt_s[melt_s$Time == i,]
      
      if (graph == 'plotly'){
        graf[[i]] <- dplyr::`%>%`(plotly::plot_ly(
          x = as.numeric(as.character(melt_s_2$X)), #melt_s_2$X,
          y = as.numeric(as.character(melt_s_2$Y)), #melt_s_2$Y,
          z = melt_s_2$Value,
          type = "heatmap",
          colorbar = list(title = "Prediction"),
          reversescale = T,
          zmin = zmin,
          zmax = zmax
        ),
        plotly::add_trace(x = KS$SFD[[1]]$coords[,1], 
                  y = KS$SFD[[1]]$coords[,2], text=sprintf("<b>%s</b>", rownames(KS$SFD[[1]]$coords)), 
                  type = 'scatter', mode = 'markers+text',
                  textposition = "center",
                  textfont = list(color = 'white',size=25,family = "Century Gothic",face="bold"),
                  marker = list(size = 50, color = 'black',opacity = 0.5,
                                line = list(
                                  color = 'rgb(231, 99, 250)',
                                  width = 8
                                )))
        )
      }
      if (graph == 'gg'){
        graf[[i]] <- ggplot2::ggplot(data = NULL,
                                     ggplot2::aes(x = as.numeric(as.character(melt_s_2$X)), #melt_s_2$X,
                                                  y = as.numeric(as.character(melt_s_2$Y))))+ #melt_s_2$Y,
          ggplot2::geom_tile(ggplot2::aes(fill = melt_s_2$Value))+
          ggplot2::labs(fill = "Prediction",title = paste("Prediction - Time = ", times[i]),
                        x = '',y = '',color = NULL,lwd = NULL)+
          ggplot2::scale_fill_viridis_c(direction = -1,limits = c(zmin,zmax)) +
          ggplot2::coord_fixed() +
          ggplot2::theme(plot.background = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         panel.border = ggplot2::element_blank(),
                         axis.line = ggplot2::element_blank(),
                         axis.text = ggplot2::element_blank())
      }
      
    }
    
    return(graf)
    
  }