# HPRD Normal - both unknown
# R: predictive R = ap
# P: predictive P = bp/(1+bp)

HPRD_Poisson <- function(far, R, P){

  lb <- max(  floor(  R*(1-P)/P - sqrt( (1/far)* R*(1-P)/(P^2) )  ), 0  )
  ub <- ceiling(  R*(1-P)/P + sqrt( (1/far) * R * (1-P)/(P^2) )  )

  # locations of ordered probabilities
  Pi <- order( stats::dnbinom( lb:ub, size=R, prob=P ), decreasing=T ) + lb - 1
  nnn <- 1
  sumprob <- 0
  diff <- 1
  E <- c()
  stopp <- 0

  while (stopp==0) {
    sumprob <- sumprob + stats::dnbinom( Pi[nnn], size=R, prob=P )
    if ( abs(sumprob - (1-far)) < diff ) {
      E <- c( E, Pi[nnn] )
      diff <- abs( sumprob - (1-far) )
      nnn <- nnn+1
    } else { stopp <- 1 }
  }


  hprd <- c( min(E), max(E) )

  return( hprd )
}



# Function for running PCC
PCC_Poisson <- function( data = NULL, s = NULL, historical_data = NULL,
                         historical_s = NULL, c = 1/2, d = 0, alpha_0 = NULL,
                         ARL_0 = 370.4, FAP = NULL, FIR = FALSE, fFIR = .99,
                         aFIR=1/8, summary_list = TRUE, PCC_PLOT = TRUE,
                         historical_data_PLOT = FALSE, pdf_report = FALSE,
                         xlab = "Observations", ylab = "Values", main="PCC Poisson")
{
  ### Initial checks before procceeding to the main body of function
  ### Mainly this chunk of code will correspond to invalid general input before running stuff
  # 'data' (i) not defined (ii) not in vector (iii) contain non-numeric value
  if ( is.null(data) ) {
    stop("'data' have not been defined")
  } else { if ( any(!is.numeric((unlist(data)))) | any(unlist(data)<0) | any(unlist(data)%%1!=0) ) stop("Invalid 'data' input")
    if ( !is.vector(data) ) stop("'data' must be in vector form")
  }
  # 'historical_data' (i) not in vector (ii) contain non-numeric value
  if ( !is.null(historical_data) ) {
    if ( any(!is.numeric((unlist(historical_data)))) | any(unlist(historical_data)<0) | any(unlist(historical_data)%%1!=0) ) stop("Invalid 'historical_data' input")
    if ( !is.vector(data) ) stop("'historical data' must be in vector form")
  }


  # 'ARL_0' (i) non-numeric (ii) negative
  if( !missing(ARL_0) ) {
    if ( length(unlist(ARL_0))>1 ) { message("More than one value for 'ARL_0', the first one will only be used")
      if ( !is.numeric(ARL_0[1]) | ARL_0<=0 ) { stop("Invalid 'ARL_0' value") } else { ARL_0 <- ARL_0[1] }
    } else { if ( !is.numeric(ARL_0) | ARL_0<=0 ) { stop("Invalid 'ARL_0' value") } }
  }

  # 'FAP' (i) non-numeric (ii) negative
  if (!missing(FAP)){
    if ( length(unlist(FAP))>1 ) { message("More than one value for 'FAP', the first one will only be used")
      if ( !is.numeric(FAP[1]) | FAP<=0 | FAP>=1 ) { stop("Invalid 'FAP' value") } else { FAP <- FAP[1] }
    } else { if ( !is.numeric(FAP) | FAP<=0 | FAP>=1 ) { stop("Invalid 'FAP' value") } }
  }

  # 'FIR' (i) logical (ii) fFIR - aFIR conditions
  if ( length(unlist(FIR))>1 ) {
    message("More than one value for 'FIR', the first one will only be used")
    if ( !is.logical(FIR[1]) ) { stop("Invalid 'FIR' value ; 'FIR' must be logical") } else { FIR <- FIR[1] }
  } else {
    if ( !is.logical(FIR) ) { stop("Invalid 'FIR' value ; 'FIR' must be logical") }
  }

  # fFIR - aFIR conditions if  FIR
  if ( FIR ) {
    if ( !missing(fFIR) ) {
      if ( length(unlist(fFIR))>1 ) {
        message("More than one value for 'fFIR', the first one will only be used")
        if ( !is.numeric(fFIR[1]) | fFIR[1]<=0 | fFIR[1]>=1 ) {
          stop("Invalid 'fFIR' value")
        } else { fFIR <- fFIR[1] }
      } else {
        if ( !is.numeric(fFIR) | fFIR<=0 | fFIR>=1 ) {
          stop("Invalid 'fFIR' value")
        }
      }
    }

    if ( !missing(aFIR) ) {
      if ( length(unlist(aFIR))>1 ) {
        message("More than one value for 'aFIR', the first one will only be used")
        if ( !is.numeric(aFIR[1]) | aFIR[1]<=0 ) {
          stop("Invalid 'aFIR' value")
        } else { aFIR <- aFIR[1] }
      } else {
        if ( !is.numeric(aFIR) | aFIR<=0 ) {
          stop("Invalid 'aFIR' value")
        }
      }
    }
  }

  ### Setting the False Alarm Probability & False Alarm Rate based on the Sidak correction
  # data length
  N <- length(data)
  # If both ARL_0 and FAP chosen
  if ( !is.null(ARL_0) & !is.null(FAP) ) {
    message("Both ARL_0 and FAP are defined as input, so ARL_0 is used by default. \nIn order to use FAP instead, set ARL_0 = NULL")
    FAR <- 1/ARL_0
    # If only FAP is chosen
  } else if ( is.null(ARL_0) & !is.null(FAP) ) {
    FAR <- 1-(1-FAP)^(1/(N-1))
    # If only ARL0 is chosen
  } else if ( !is.null(ARL_0) & is.null(FAP) ){
    FAR <- 1/ARL_0
  }
  # If FIR PCC is chosen - default value for f=0.99
  if ( FIR ) {
    tf <- 1:N
    Afir <- c((  1- (1-fFIR)^(1+aFIR*(tf-1)) ) )
    FAR <- 1-(1-FAR)*Afir
  }


  ###############################################################
  ###############################################################
  ## START (1) Only this bit changes from function to function ##
  ## In some cases 'data' and 'historical_data' restrictions     ##
  ## change as well at the beginning of the function           ##
  ###############################################################
  ###############################################################

  # Prior parameter input (i) more than one value for parameters (ii) non-numeric input
  if ( is.null(s) ) { s <- rep(1, times=length(data))
  } else {
    if ( !is.vector(s) ) { stop("'rates - s' must be in vector form")
    } else {
      if ( length(s)!=length(data) ) { stop("Vector of 'rates - s' must have the same length as 'data'")
      } else {
        if ( any(!is.numeric((unlist(s)))) ) { stop("Invalid 'rate - s' input")
        } else { if( any((unlist(s)<=0)) ) stop("Invalid 'rate - s' input, s must be positive") }
      }
    }
  }

  if( !missing(c) ) {
  if ( length(unlist(c))>1 ) { message("More than one value for 'c', the first one will only be used")
    if ( !is.numeric(c) | c<0 ) { stop("Invalid 'c' value") } else { c <- c[1] }
  } else { if ( !is.numeric(c) | c<0 ) { stop("Invalid 'c' value") } }
  }

  if( !missing(d) ) {
  if ( length(unlist(d))>1 ) { message("More than one value for 'd', the first one will only be used")
    if ( !is.numeric(d) | d<0 ) { stop("Invalid 'd' value") } else { d <- d[1] }
  } else { if ( !is.numeric(d) | d<0 ) { stop("Invalid 'd' value") } }
  }

  ### Main body of function - PCC illustration - USING FAR (or FAP equivelantly)
  ## Histotic data and processing
  if ( !is.null(historical_data) ){

    if ( is.null(historical_s) ) { historical_s <- rep(1, times=length(data))
    } else {
      if ( !is.vector(historical_s) ) { stop("'historical rates - historical_s' must be in vector form")
      } else {
        if ( length(historical_s)!=length(data) ) { stop("Vector of 'historical_s rates - s' must have the same length as 'historical_data'")
        } else {
          if ( any(!is.numeric((unlist(historical_s)))) ) { stop("Invalid 'historical rates - historical_s' input")
          } else { if( any((unlist(historical_s)<=0)) ) stop("Invalid 'historical rates - historical_s' input, s must be positive") }
        }
      }
    }

  N_historicaldata <- length(historical_data)
    # Check about alpha_0
    # If no chosen value for alpha_0 use default setting
    if (is.null(alpha_0)) { alpha_0 <-1/N_historicaldata
    } else {
      if ( length(unlist(alpha_0))>1 ) {
        message("More than one value for 'alpha_0', the first one will only be used")
        if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1) { stop("Invalid 'alpha_0' value")
        } else { alpha_0 <- alpha_0[1] }
      } else { if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1 ) { stop("Invalid 'alpha_0' value") } }
    }
    # Process historical data
    # Power Prior parameters
    c_PowerP  <- c + alpha_0*sum(historical_data)
    d_PowerP  <- d + alpha_0*sum(s)
    # Keep similar notation as input
    c <- c_PowerP ; d <- d_PowerP
  }

  ### PCC implementation
  # Sum of observations
  dataSum <- cumsum(data)[seq(1,length(data)-1)]
  # Sum of rates
  dataRates <- cumsum(s)[seq(1,length(data)-1)]
  # Predictive distribution parameters
  R_Pred <- c + dataSum
  P_Pred <- ( d + dataRates ) / ( d + dataRates + s[seq(2, length(data))] )

  # Control limits
  if (!FIR) { edges <- t( mapply( function(Rind, Pind, FARind=FAR) {
                                    do.call(HPRD_Poisson, list(far=FARind, R=Rind, P=Pind))
                                  }, Rind=R_Pred, Pind=P_Pred) )
  } else { edges <- t( mapply( function(Rind, Pind, FARind) {
                                do.call(HPRD_Poisson, list(far=FARind, R=Rind, P=Pind))
                               }, Rind=R_Pred, Pind=P_Pred, FARind=FAR) )
  }

  edges <- rbind( c(NA, NA), edges )


  ####################################################################
  ####################################################################
  ## END (1) Only the above bit changes from function to function   ##
  ####################################################################
  ####################################################################


  ## Output
  { # Construction of 'In' and 'Out' of control column for return results
    States <- rep("", times=N)
    States[ifelse(data < edges[, 1], TRUE, FALSE)] <- "Alarm (LL)" ; States[ifelse(data > edges[, 2], TRUE, FALSE)] <- "Alarm (UL)"
    # Return results
    PCC_summary <- data.frame(  data=data, HPrD_LL=edges[, 1], HPrD_UL=edges[, 2], Alarms=States )   }

  ## Dynamic recalculation of PCC plot's y axis
  # PCC y axis limits allowance
  Ratio <- (PCC_summary$HPrD_UL-PCC_summary$HPrD_LL)/min(PCC_summary$HPrD_UL-PCC_summary$HPrD_LL, na.rm = T)

  # Y axis limits
  AdjustedYlim <- c(min(PCC_summary$data, PCC_summary$HPrD_LL[which(Ratio<=2.5)], na.rm=T),
                    max(PCC_summary$data, PCC_summary$HPrD_UL[which(Ratio<=2.5)], na.rm=T))


  ### Output of function
  ## PCC plot
  if ( PCC_PLOT ) {
    # Creation of PCC plot
    PCC_PlotSummary <- cbind(Indices=1:N, PCC_summary)
    PCC <- ggplot2::ggplot(PCC_PlotSummary, ggplot2::aes(Indices, data)) +
      ggplot2::geom_line(ggplot2::aes(x=Indices, y=data), na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(x=Indices, y=HPrD_UL), color="red", linetype="solid", size=1, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(x=Indices, y=HPrD_LL), color="red", linetype="solid", size=1, na.rm = TRUE) +
      ggplot2::geom_ribbon(ggplot2::aes(x=Indices, ymin=HPrD_UL, ymax=HPrD_LL, fill=TRUE), alpha=0.25, show.legend=FALSE) +
      ggplot2::scale_fill_manual(values=c("TRUE"="green")) +
      ggplot2::geom_point(ggplot2::aes(group=Indices, color=as.factor(Alarms), stroke = 1.5), show.legend=FALSE, na.rm = TRUE) +
      ggplot2::scale_color_manual(values=c("black", "red", "red"), na.value = "black") +
      ggplot2::coord_cartesian(ylim = AdjustedYlim) +
      ggplot2::labs(title = main, x = xlab, y = ylab) +
      ggplot2::theme(legend.position = "top",
                     legend.title = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(colour = "black", size = 0.5, linetype = "solid"),
                     panel.background = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5)
      )
    # Creation of PCC plot if historical data are chosen to be on the plot
    if ( !is.null(historical_data) & historical_data_PLOT ) {
      PCC_summary_historicaldata <- data.frame(  data=c(historical_data, data), HPrD_LL=c(rep(NA, times=N_historicaldata), edges[, 1]),
                                               HPrD_UL=c(rep(NA, times=N_historicaldata), edges[, 2]), Alarms=c(rep("", times=N_historicaldata), States) )
      PCC_PlotSummary <- cbind(Indices=c(-N_historicaldata:(-1), 1:N), TypeOfdata=c(rep("Historical", times=N_historicaldata), rep("Current", times=N)), PCC_summary_historicaldata )
      PCC_historical <- ggplot2::ggplot(PCC_PlotSummary, ggplot2::aes(Indices, data)) +
        ggplot2::geom_line(ggplot2::aes(x=Indices, y=data, linetype = as.factor(TypeOfdata)), na.rm = TRUE) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, y = min(HPrD_LL, na.rm=TRUE), xend = 0, yend = max(HPrD_UL, na.rm=TRUE))) +
        ggplot2::geom_line(ggplot2::aes(x=Indices, y=HPrD_UL), color="red", linetype="solid", size=1, na.rm = TRUE) +
        ggplot2::geom_line(ggplot2::aes(x=Indices, y=HPrD_LL), color="red", linetype="solid", size=1, na.rm = TRUE) +
        ggplot2::geom_ribbon(ggplot2::aes(x=Indices, ymin=HPrD_UL, ymax=HPrD_LL, fill=TRUE), alpha=0.25, show.legend=FALSE) +
        ggplot2::scale_fill_manual(values=c("TRUE"="green")) +
        ggplot2::geom_point(ggplot2::aes(group=Indices, shape=as.factor(TypeOfdata), color=as.factor(Alarms), stroke = 1.5), show.legend=FALSE, na.rm = TRUE) +
        ggplot2::scale_color_manual(values=c("black", "red", "red"), na.value = "black") +
        ggplot2::scale_linetype_manual(values=c("Historical"="dotted", "Current"="solid")) +
        ggplot2::scale_shape_manual(values=c("Historical"=1, "Current"=19)) +
        ggplot2::coord_cartesian(ylim = AdjustedYlim) +
        ggplot2::labs(title = main, x = xlab, y = ylab) +
        ggplot2::theme(legend.position = "top",
                       legend.title = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black", size = 0.5, linetype = "solid"),
                       panel.background = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)
        )
      print(PCC_historical)

    } else { print(PCC) }

  }
  # List of results
  if ( summary_list ) { print(PCC_summary) }

  # List of results return in pdf
  if ( pdf_report ) {

    # save pdf
    grDevices::pdf(
      paste0( "PCC_results_", paste0( unlist(strsplit(date(), " "))[c(1,2,3,5)], collapse = "_" ), "_",
              paste0( unlist(strsplit( unlist(strsplit(date(), " "))[4], ":" )), collapse = "." ),
              ".pdf" ),
      height = 8.264, width = 11.694)

    # PCC plot on pdf
    print(PCC)

    # Results matrix on pdf
    # Chunk of code to split results matrix to different pages - Set a default number based on pdf height/width
    NRowsPerPage <- 25
    if(NRowsPerPage > nrow(PCC_summary)){ FloatingRow <- nrow(PCC_summary) } else { FloatingRow <- NRowsPerPage }
    sapply(1:ceiling(nrow(PCC_summary)/NRowsPerPage), function(index) {
      if (index==1) { StartingRow <<- 1 }
      grid::grid.newpage()
      gridExtra::grid.table(PCC_summary[StartingRow:FloatingRow, ])
      StartingRow <<- FloatingRow + 1
      if( sum(NRowsPerPage, FloatingRow) < nrow(PCC_summary)){ FloatingRow <<-  NRowsPerPage + FloatingRow } else { FloatingRow <<- nrow(PCC_summary) }
    })

    grDevices::dev.off()

  }


}






