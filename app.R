###  ---------
## name: mvn_sim.R
## Author: Kyle N. Payne
## Purpose: simulate sets of multivariate normal random variables in
## order to test the effects of autocorrelated-ness in the 
## random effects structure data.
###  ---------


source("helpers.R")
require(lme4)
require(plyr)
require(shiny)
require(fdrtool)
require(mvtnorm)
require(reshape2)
require(eegkit)
data(eegcoord)


### ---- list of inputs and outputs for the program
### ---- inputs:
### ---- min_sill
### ---- max_sill
### ---- min_decay
### ---- max_decay
### ---- num_steps (decay and sill evaluation points)
### ---- region (region of the brain)
### ---- effect (effect size)
### ---- nugget
### ---- alpha 
### ---- num_subs
### ---- num_iter
### ---- samp_size
### ---- outputs: 
### ---- dot_map

### ---- server part of the program.
server <- function(input, output){

 output$dot_map <- renderPlot({
    ## make sure that computation starts with goButton
    input$goButton
   
    min_sill <- input$min_sill ; max_sill <- input$max_sill
    min_decay <- input$min_decay ; max_decay <- input$max_decay
    num_steps <- input$num_steps
    num_iter <- input$num_iter
    region <- input$region
    effect <- input$effect
    alpha <- input$alpha
    D <- input$num_subs
    N <- input$samp_size
    q <- nrow(eegcoord)
    nug <- input$nugget
    silly <- seq(from=min_sill, to=max_sill, length.out=num_steps)
    decay <- seq(from=min_decay, to=max_decay, length.out=num_steps)
    locs <- c("xproj","yproj")
    
    ## ---- create a matrix of lists that contains
    ## ---- all of the different values of the psill and
    ## ---- of the decay
    
    params <- matrix(rep(list(), length(silly)*length(decay)), nrow=length(silly), ncol=length(decay))
    p_values <- matrix(rep(list(), length(silly)*length(decay)), nrow=length(silly), ncol=length(decay))
    
    ## fill the paramter matrix of lists
    for(s in 1:length(silly)){
      for(d in 1:length(decay)){
        params[s,d][[1]] <- c(silly[s], decay[d])
      }
    }  
    
    withProgress(message = 'Making plot', value = 0, {
    for(s in 1:length(silly)){
      for(d in 1:length(decay)){
        incProgress(1/num_steps^2, detail = "Sorry, this takes awhile")
          if(input$covar_type=="expo"){
            ### ---- returns the exponential covariance
            covar <- cov_exp(nugget = nug, 
                             psill = silly[s], 
                             decay= decay[d],
                             locations = locs,
                             data = eegcoord
            )
          }
          if(input$covar_type=="sphere"){
            ### ---- returns the spherical covariance
            covar <- cov_sphere(nugget = nug,
                                psill = silly[s],
                                decay = decay[d],
                                locations = locs,
                                data=eegcoord)
          }
          if(input$covar_type=="gauss"){
            ### ---- returns the gaussian covariance
            covar <- cov_gauss(nugget = nug,
                               psill=silly[s],
                               decay=decay[d],
                               locations=locs,
                               data=eegcoord)
          }
          
          ## ----- check to make sure that the 
          ## ---- covariance matrix is positive definite
          ## ---- if not, set the p-value list to 
          ## ---- NA and move onto the next iteration
          if(any(is.na(covar))){
            p_values <- NA
            next
          }
          ## ----- empty data frame to hold the simulations for
          ## ----- all numbers of participants
          data_sim_final <- data.frame(chan=character(),
                                       c2=numeric(),
                                       mean=numeric(),
                                       xproj=numeric(),
                                       yproj=numeric(),
                                       subject=numeric())
          
        ## ----- we use this dataset for coordinates of eegcap locations
      eegcoord_edit <- cbind(eegcoord, row.names(eegcoord))
      eegcoord_edit <- eegcoord_edit[,-c(1,2,3)]
      colnames(eegcoord_edit)[ncol(eegcoord_edit)] <- "chan"
      ## ---- initialize a NULL vector that we will fill 
      ## ---- with the proportion of times that the
      ## ---- corresponding pval from the t-test
      ## ---- of the c2 effect is less than 0.05
      p_vals <- NULL
      
      for(iter in 1:num_iter){   
        ## for each iteration, calculate the 
        ## proportion of times that the 
        ## p_value is non-significant
        for(subj in 1:D){
          ## create means vectors 
          mu_1 <- rnorm(q, 0 ,1) 
          ## mu_1 <- rep.row(mu_1, N)
          
          mu_2 <- rnorm(q, 0, 1)
          ## mu_2 <- rep.row(mu_2, N)
          
          ### ---- induce a pocket of electrical activity
          ### ---- that is clustered about the selected electrode sites
          electrodes <- row.names(eegcoord)
          elect_labels <- which(substr(electrodes, 1,1) == region)
          mu_2[elect_labels] <- rnorm(length(elect_labels), -1*effect, 2)
          ## ---- selecting a variance of 2 is a completely arbitrary choice
          ## ---- and I'll need to change it in subsequent versions of the program
      
          ## ---- create fake EEG data with different mean shifts
          ## ---- but with identical covariances
          y_1 <- rmvnorm(N, mean = mu_1, sigma=covar) 
          y_2 <- rmvnorm(N, mean = mu_2, sigma=covar) 
          
          ## ---- fake meanAmpEEG electrode data
          data_sim <- rbind(y_1, y_2)
          data_sim <- cbind(data_sim, c(rep(0, N), rep(1, N)))
          data_sim_df <- data.frame(data_sim)
          colnames(data_sim_df) <- c(row.names(eegcoord), "c2")
          data_sim_df_m <- melt(data_sim_df, id=.(c2))
          colnames(data_sim_df_m) <- c("c2","chan", "meanAmp")
          data_sim_df_m_merg <- merge(data_sim_df_m, eegcoord_edit, by="chan")
          data_sim_df_m_merg <- cbind(data_sim_df_m_merg, subject=subj)
          data_sim_final <- rbind(data_sim_final, data_sim_df_m_merg) 
          rm(data_sim_df_m_merg)
        }  
    
      model_mixed <- lmer(meanAmp ~ c2 + (1|chan) + (1|subject), 
                          data=data_sim_final)
      
      ### get the t-statistic values from the summary function
      
      t_stats <- summary(model_mixed)$coefficients[,3]
      
      ## normal approximation appears to be slightly anti-conservative.
      p_stats <- 1-pnorm(abs(t_stats))
      
      c2_p <- p_stats['C2']
      p_vals <- c(p_vals, c2_p)
      }
      ## ---- put the p_values into a list, within a matrix 
      ## ---- whose entries are determined by the 
      ## ---- number of sill and decay parameters
      p_values[s,d][[1]] <- p_vals
      }
    }
    }) ## end of with progress loop
    ## ---- outside of sill and decay loop.
    ## --- let's look at the proportion of 
    ## --- signif p_vals
    prop <- matrix(NA, nrow=length(silly), ncol=length(decay))
    for(s in 1:length(silly)){
      for(d in 1:length(decay)){
          prop[s,d][[1]] <- mean(p_values[s,d] >= alpha)
        }
    }
    colnames(prop) <- silly
    row.names(prop) <- decay
    prop_melt <- melt(prop)
    colnames(prop_melt) <- c("psill", "decay", "prop_sig_pvals")
    
    ## ---- dot_map will be a ggplot where the color of the dot
    ## ---- describes the proportion of significant p_values 
    dot_map <- ggplot(aes(x=decay, y=psill), data=prop_melt) +
               geom_point(aes(color=prop_sig_pvals))
    dot_map
    })
}
### ---- end of server part
###
###
###
###
### ----- user interface part
### ---- list of inputs and outputs for the program
### ---- inputs:
### ---- min_sill
### ---- max_sill
### ---- min_decay
### ---- max_decay
### ---- num_steps (decay and sill evaluation points)
### ---- region (region of the brain)
### ---- effect (effect size)
### ---- nugget
### ---- alpha 
### ---- outputs: 
### ---- dot_map
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      numericInput("min_sill", label="Min psill",
                   min=0.01, max=100, value=10, step=0.01),      
      numericInput("max_sill", label="Max psill",
                   min=0.01, max=100, value=10, step=0.01),
      numericInput("min_decay", label="Min Decay",
                   min=0.01, max=100, value=10, step=0.01),
      numericInput("max_decay", label="Max Decay",
                   min=0.01, max=100, value=10, step=1),
      numericInput("num_steps", label="Number of Steps",
                   min=1, max=1000, value=10, step=1), 
      numericInput("num_iter", label="Number of Iterations",
                   min=1, max=1000, value=10, step=1),
      numericInput("num_subs", label="Number of Subjects",
                   min=1, max=100, value=10, step=1),
      numericInput("samp_size", label = "Number of Trials",
                   min = 0.001, max = 10000, value = 1, step = 1),
      numericInput("nugget", label = "Nugget",
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
                    "occipital" = "O")),
      actionButton("goButton", "Start")),
    mainPanel(
      plotOutput("dot_map")
    )
  )
)

shinyApp(ui = ui, server = server)
