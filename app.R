#### -----
## app.R 
## Kyle N. Payne
#### ----

library(shiny)
source("helpers.R")

library(mvtnorm)
library(reshape2)
library(plyr)
library(eegkit)
library(lme4)
data(eegcoord)

## ---- Server 

server <- function(input, output){
  output$scalpPlot <- renderPlot({
      ### ------------------------------------------------------------ ###
      ### --------- Generating spatially correlated random fields ---- ###
      ### ------------------------------------------------------------ ###
      D <- input$num_subs
      N <- input$samp_size
      q <- nrow(eegcoord)
      nug <- input$nugget
      silly <- input$part_sill
      phi <- input$phi
      locs <- c("xproj","yproj")
      
      if(input$covar_type=="expo"){
        ### ---- returns the exponential covariance
        covar <- cov_exp(nugget = nug, 
                             psill = silly, 
                             decay= phi,
                             locations = locs,
                             data = eegcoord
        )
      }
      if(input$covar_type=="sphere"){
        covar <- cov_sphere(nugget)
      }
      if(input$covar_type=="gauss"){
        covar <- cov_mat()
      }
      ## ----- empty data frame to hold the simulations for
      ## ----- all numbers of participants
      data_sim_final <- data.frame(chan=character(),
                                  c2=numeric(),
                                  mean=numeric(),
                                  xproj=numeric(),
                                  yproj=numeric(),
                                  subject=numeric())
      ## ----- we use this dataset to 
      eegcoord_edit <- cbind(eegcoord, row.names(eegcoord))
      eegcoord_edit <- eegcoord_edit[,-c(1,2,3)]
      colnames(eegcoord_edit)[ncol(eegcoord_edit)] <- "chan"
      for(subj in 1:D){
        ## create means vectors 
        mu_1 <- rnorm(q, 0 ,1) 
        ## mu_1 <- rep.row(mu_1, N)
        
        mu_2 <- rnorm(q, 0, 1)
        ## mu_2 <- rep.row(mu_2, N)
        
        ### ---- induce a pocket of electrical activity
        ### ---- that is clustered about the centroParietal sites
        electrodes <- row.names(eegcoord)
        elect_labels <- which(substr(electrodes, 1,1) == input$region)
        mu_2[elect_labels] <- rnorm(length(elect_labels), -1*input$effect, 2)
        
        ## ---- create fake EEG data with different mean shifts
        ## ---- but with identical covariances
        y_1 <- rmvnorm(N, mean = mu_1, sigma=covar) 
        y_2 <- rmvnorm(N, mean = mu_2, sigma=covar) 
        
        ## ---- fake meanAmpEEG electrode data
        data_sim <- rbind(y_1, y_2)
        data_sim <- cbind(data_sim, c(rep(0, N), rep(1, N)))
        data_sim_df <- data.frame(data_sim)
        colnames(data_sim_df) <- c(row.names(eegcoord), "c2")
        data_sim_df_m <- melt(data_sim_df, .(c2))
        colnames(data_sim_df_m) <- c("c2","chan", "meanAmp")
        data_sim_df_m_merg <- merge(data_sim_df_m, eegcoord_edit, by="chan")
        data_sim_df_m_merg <- cbind(data_sim_df_m_merg, subject=subj)
        ## ---- let's look at a scalp-map lattice of the meanAmpEEG
        ## ---- averages for each channel and across the two
        ## ---- conditions

        data_sim_final <- rbind(data_sim_final, data_sim_df_m_merg) 
    }  
    data_sim_means <- ddply(data_sim_final,
                          .(c2, chan, subject), summarize,
                          mean=mean(meanAmp),
                          xproj = max(xproj),
                          yproj = max(yproj))
    

    
    for(j in 1:D){
      for(i in 1:2){
        ## ---- for each participant, and each condition
        eegspace(
        ## ---- grab the x-y projected location  
        data_sim_means[which(data_sim_means$c2 == i-1
                  & data_sim_means$subject == j),
                                c("xproj", "yproj")],
        ## ---- grab the means        
        voltage = data_sim_means[which(data_sim_means$c2 == i-1
                  & data_sim_means$subject == j),
                                "mean"], 
        vlim=c(min(data_sim_means$mean),max(data_sim_means$mean)),
        plotaxes=FALSE,
        cex.axis=.5,
        colorbar=FALSE)
      }
    }
  })
}
## ---- User Interface

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
    numericInput("num_subs", label="Number of Subjects",
                 min=1, max=100, value=10, step=1),
    numericInput("samp_size", label = "Number of Trials",
                 min = 0.001, max = 10000, value = 1, step = 1),
    numericInput("part_sill", label = "Partial Sill",
                 min = 0.001, max = 500, value = 1, step = 0.01),
    numericInput("nugget", label = "Nugget",
                 min = 0.001, max = 500, value = 1, step = 0.01),
    numericInput("phi", label = "Decay",
                 min = 0.001, max = 500, value = 1, step = 0.01),
    sliderInput("effect", label = "Effect Size",
                min = 0.001, max = 500, value = 1, step = 0.01),
    selectInput("covar_type", "Type of Spatial Covariance",
                c("Exponential" = "expo",
                  "Spherical" = "sphere",
                  "Gaussian" = "gauss")),
    selectInput("region", "Electrode Regions",
                c("frontal" = "F", 
                  "anterior" = "A", 
                  "temporal" = "T", 
                  "central" = "C", 
                  "parietal" = "P", 
                  "occipital" = "O"))),
    mainPanel(
      plotOutput("scalpPlot")
      )
  )
)

shinyApp(ui = ui, server = server)