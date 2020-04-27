# HPRD Normal - both unknown
HPRD_Normal_BothUnknown <- function(far, Mpr, Lp, Ap, Bp){

  hprd <- c( Mpr + stats::qt( far/2, df=2*ifelse(Ap > 0, Ap, NA) ) * sqrt( ((Bp*(Lp+1)) / (Ap*Lp)) ),
             Mpr + stats::qt( 1-far/2, df=2*ifelse(Ap > 0, Ap, NA) ) * sqrt( ((Bp*(Lp+1)) / (Ap*Lp)) )
            )
  return( hprd )

}


# Function for running PCC
PCC_Normal_BothUnknown <- function( data = NULL, historical_data = NULL,
                                    mu0 = 0, l0 = 0, a0 = -1/2, b0 = 0, alpha_0 = NULL,
                                    ARL_0 = 370.4, FAP = NULL, FIR = FALSE, fFIR = .99,
                                    aFIR=1/8, summary_list = TRUE, PCC_PLOT = TRUE,
                                    historical_data_PLOT = FALSE, pdf_report = FALSE,
                                    xlab = "Observations", ylab = "Values", main="PCC Normal")
{
  ### Initial checks before procceeding to the main body of function
  ### Mainly this chunk of code will correspond to invalid general input before running stuff
  # 'data' (i) not defined (ii) not in vector (iii) contain non-numeric value
  if ( is.null(data) ) {
    stop("'data' have not been defined")
  } else { if ( any(!is.numeric((unlist(data)))) ) stop("Invalid 'data' input")
    if ( !is.vector(data) ) stop("'data' must be in vector form")
  }
  # 'historical_data' (i) not in vector (ii) contain non-numeric value
  if ( !is.null(historical_data) ) {
    if ( any(!is.numeric((unlist(historical_data)))) ) stop("Invalid 'historical_data' input")
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
  ## In some cases 'data' and 'historical_data' restrictions   ##
  ## change as well at the beginning of the function           ##
  ###############################################################
  ###############################################################

  # Prior parameter input (i) more than one value for parameters (ii) non-numeric input
  if( !missing(mu0) ) {
  if ( length(unlist(mu0))>1 ) { message("More than one value for 'mu0', the first one will only be used")
    if ( !is.numeric(mu0[1]) ) { stop("Invalid 'mu0' value") } else { mu0 <- mu0[1] }
  } else { if ( !is.numeric(mu0) ) { stop("Invalid 'mu0' value") } }
  }

  if( !missing(l0) ) {
  if ( length(unlist(l0))>1 ) { message("More than one value for 'l0', the first one will only be used")
    if ( !is.numeric(l0) | l0<0 ) { stop("Invalid 'l0' value") } else { l0 <- l0[1] }
  } else { if ( !is.numeric(l0) | l0<0 ) { stop("Invalid 'l0' value") } }
  }

  if( !missing(a0) ) {
  if ( length(unlist(a0))>1 ) { message("More than one value for 'a0', the first one will only be used")
    if ( !is.numeric(a0) | a0<0 ) { stop("Invalid 'a0' value") } else { a0 <- a0[1] }
  } else { if ( !is.numeric(a0) | a0<0 ) { stop("Invalid 'a0' value") } }
  }

  if( !missing(b0) ) {
  if ( length(unlist(b0))>1 ) { message("More than one value for 'b0', the first one will only be used")
    if ( !is.numeric(b0) | b0<0 ) { stop("Invalid 'b0' value") } else { b0 <- b0[1] }
  } else { if ( !is.numeric(b0) | b0<0 ) { stop("Invalid 'b0' value") } }
  }

  ### Main body of function - PCC illustration - USING FAR (or FAP equivelantly)
  ## Histotic data and processing
  if ( !is.null(historical_data) ){
    N_historicaldata <- length(historical_data)
    # Check about alpha_0
    # If no chosen value for alpha_0 use default setting
    if (is.null(alpha_0)) { alpha_0 <-1/N_historicaldata
    } else {
      if ( length(unlist(alpha_0))>1 ) {
        message("More than one value for 'alpha_0', the first one will only be used")
        if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1) { stop("Invalid 'alpha_0' value")
        } else { if ( !is.numeric(alpha_0) | alpha_0<0 | alpha_0>1 ) { stop("Invalid 'alpha_0' value") } }
      }
    }
    # Process historical data
    # Power Prior parameters
    mu0_PowerP <- (l0*mu0+alpha_0*sum(historical_data))/(l0+alpha_0*N_historicaldata)
    l0_PowerP  <- l0 + alpha_0*N_historicaldata
    a0_PowerP  <- a0 + alpha_0*N_historicaldata/2
    b0_PowerP  <- b0 + alpha_0*sum(historical_data^2)/2 + (l0*mu0^2)/2-((alpha_0*sum(historical_data)+l0*mu0)^2)/(2*(l0+alpha_0*N_historicaldata))
    # Keep similar notation as input
    mu0 <- mu0_PowerP ; l0 <- l0_PowerP ; a0 <- a0_PowerP ; b0 <- b0_PowerP
  }

  ### PCC implementation
  # Sum of observations
  dataSum <- cumsum(data)[seq(1,length(data))]
  # Sum of squared observations
  dataSumSquared <- cumsum(data^2)[seq(1,length(data))]
  # Predictive distribution parameters
  mu0_Pred <- (l0*mu0+dataSum)/(l0+(1:N))
  l0_Pred <- l0+(1:N)
  a0_Pred <- a0+(1:N)/2
  b0_Pred <- b0+dataSumSquared/2+(l0*mu0^2)/2-((dataSum+l0*mu0)^2/(2*(l0+(1:N))))


  # Control limits
  if (!FIR) { edges <- t( mapply( function(MU0, L0, A0, B0, FD=FAR) { do.call(HPRD_Normal_BothUnknown, list(far=FD, Mpr=MU0, Lp=L0, Ap=A0, Bp=B0)) }
                                  , MU0=mu0_Pred, L0=l0_Pred, A0=a0_Pred, B0=b0_Pred) )
  } else { edges <- t( mapply( function(MU0, L0, A0, B0, FD) { do.call(HPRD_Normal_BothUnknown, list(far=FD, Mpr=MU0, Lp=L0, Ap=A0, Bp=B0)) }
                               , MU0=mu0_Pred, L0=l0_Pred, A0=a0_Pred, B0=b0_Pred, FD=FAR) )
  }

  edges <- rbind( c(NA, NA), edges )
  edges <- edges[-nrow(edges), ]


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






