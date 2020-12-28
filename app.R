library(shiny)
ui<- fluidPage(
  #alles was auf dem user interface zu sehen sein soll
  #all parts of UI seperated by commas!
  #pic start size: 826*64
  #tag$ calls html command tags
  tags$img(height=85.18, width=1115, src="TEE_logo_header.png",
           style="display: block; margin-left: auto; margin-right: auto;"), #puts picture
  tags$h1("The real pain of climate change:",style="color:#063d79;", 
          style="text-align:center;", style="font-weight:bold;"),  #first header, giving colour in hexcode
  h3("It is increasing variance in environmental conditions that makes natural populations suffer most", style="text-align:center"), #third header
  tags$hr(), #makes horizontal line
  h6("Charlotte S. Sieger, Marleen M.P. Cobben, Thomas Hovestadt", 
     style="text-align:left"), #sixth header, centered
  tags$hr(), #horizontal line
  #creating two columns, with width 6 (out of 12)
  column(6, wellPanel(#makes a panel and HTML() makes html text, tag$p makes new paragraph
    tags$p(HTML("Most individuals face variable environmental conditions. 
           Strategies to cope with such variation are e.g. bet-hedging, 
           dispersal, or tolerance. With a systematic trend in 
           temperature or other environmental characteristics, e.g. 
           under climate change, species also experience selection 
           pressure towards a changing environmental optimum. 
           Here the evolution of niche optimum h&#772 and width 
           (tolerance) g&#772 in isolated populations under different scenarios
           is simulated.
           Tolerance trades off against maximum fertility (fitness)."), 
           style="text-align:justify;"), #justified text
    tags$p(HTML("You can select the scenario by selecting an annual increase for
              either the environmental mean H&#772, standard deviation &sigma; or both; 
            as well as the trade-off strength &alpha;, to explore the population's 
            demography and evolution."),
            style="text-align:justify"), #new paragraph, justified
    tags$p("In the graphs, a green line marks 
           the start of the chosen scenario, while a grey line shows 
           the carrying capacity above the population size.",
           style="text-align:justify"), #new paragraph, justified
    tags$hr(), #horizontal line
    radioButtons(inputId = "checkalpha", choices = c("strong", "intermediate", "weak"), 
                 selected="intermediate", label = HTML("Trade-off strength &alpha;:"), inline = TRUE), #creates radio buttons that can react, but don't yet
    sliderInput(inputId = "meanslider", label=HTML("Annual increase in mean &delta;(H&#772)"),  #creates slider
                min = 0, max = 0.03, value = 0.0), 
    sliderInput(inputId = "sdslider", label=HTML("Annual increase in sd &delta;(&sigma;)"), 
                min = 0, max = 0.03, value = 0.0)
  ), #ends well panel
  actionButton("start", "Start Simulation Again") #makes a button that starts an action
  ),
  column(6, plotOutput("environment", width = 600, height = 900)) #creates space for the plot
  
)


server<- function(input, output){
  #alles was input oder output sein soll
  
  observeEvent(input$start, { 

  output$environment<-renderPlot({
    #start simulation 
    T0=0 # start mean temp
    Tch= reactive({input$meanslider}) #0.0 # annual increase in T, read in from slider 1
    sd0=1
    Tsd= reactive({input$sdslider}) # 0.01 #annual increase in sd, read in from slider 2
    Tini=500  #burn-in length
    Ttrend=300  #scenario-length
    sigmaT<-c(rep(sd0,  Tini), sd0+(1:Ttrend)*Tsd())
    Climate=c(rep(T0,  Tini), T0+(1:Ttrend)*Tch())
    ranClimate= Climate+rnorm(length( Climate),0,sigmaT)
    trem=400 # initial data to skip from presentation
    # population parameters
    K=10000
    R0=10
    h0=T0
    if(input$checkalpha=="weak"){ #checks inout from radiobuttons and sets value accordingly
      alpha=4
    }
    if(input$checkalpha=="intermediate"){
      alpha=2
    }
    if(input$checkalpha=="strong"){
      alpha=1
    }
    par.tradeoff=1/alpha
    g2=sqrt(1/par.tradeoff) #tolerance median
    sigma.g=0.1 # initial sd for tolerance trait
    #mut.h=0.03 #mutation width h; all individuals mutate slightly by adding from standar norm
    #mut.g=0.03 # tolerance par. - all mutate by multiplying from range 1-mut.t to 1+mut.t
    mut.h=0.5 # Vers. Charlotte
    mut.g=0.5
    mut.p=0.001
    a=( R0-1)/( K* R0) #survival constant in pop model
    ## initialization
    v.meanH=numeric(length( Climate))
    v.N=numeric(length( Climate)); 
    v.fert=numeric(length( Climate))
    v.tol=numeric(length( Climate))
    v.opt=numeric(length( Climate))
    
    population=rnorm( K,  h0, 1) # K individuals with random values for optimum centered around T0
    tolerance =rlnorm( K, log(g2), sigma.g) # ... and their tolerance parameter
    
    N=length(population)
    maladapt=((population- ranClimate[1])^2)/(tolerance^2)
    tradeoff=exp(-(tolerance^2*par.tradeoff^2*0.5))
    fert= R0*tradeoff*exp(-maladapt)
    v.N[1]=N
    v.meanH[1]=mean(population)
    v.fert[1]=mean(fert)
    v.tol[1]=mean(tolerance)
    v.opt[1]=mean(population)
    
    for(i in 2:length( ranClimate))
    {
      Offspring=rpois(length(population), fert) # random number of eggs for each ind.
      nOff=sum(Offspring)
      psurv=1/(1+a*nOff)
      survOff=rbinom(N, Offspring, psurv)
      population=rep(population, survOff)
      tolerance=rep(tolerance, survOff)
      if(length(population!=0)){
        N=length(population)
      } #ends if
      else{
        N=0
      } #ends else
      #  
      #mutations
      population=population+rnorm(N, 0, mut.h) # mutate h trait
      tolerance=tolerance*runif(N, 1-mut.g, 1+mut.g) # mutate g trait
      
      muts1=rbinom(N, 1, mut.p)
      population=population+runif(N, 1-mut.h, 1+mut.h)*muts1
      muts2=rbinom(N, 1, mut.p)
      tolerance=tolerance*(1+muts2*(runif(N, -mut.g, mut.g)))
      #  
      #  tolerance=sample(tolerance) # breaks linkage!
      #  
      tradeoff=exp(-(0.5*tolerance^2*par.tradeoff^2))
      maladapt=((population- ranClimate[i])^2)/(tolerance^2)
      fert= R0*tradeoff*exp(-maladapt)
      v.N[i]=N
      v.meanH[i]=mean(population)
      v.fert[i]=mean(fert)
      v.tol[i]=mean(tolerance)
      v.opt[i]=mean(population)
    } #ends loop
    #dev.new(height=10)
    #plots the results in the assigned space
    par(mfrow=c(4,1), ps=20, mar=c(5,5,1,1), mar=c(2,5,1,1), oma=c(3,1,1,1))
    plot(na.omit( ranClimate[-(1: trem)]), type="l", bty="l",
         xlab='', ylab=expression(paste('Environmental mean ', bar(H))))
    abline(v= Tini- trem, col="lightgreen", lty=5, lwd=2)
    lines( Climate[-(1: trem)], col='#063d79', lwd=2)
    # lines(v.meanH[-(1:trem)], col="red3", lwd=2)
    
    plot(na.omit(v.N[-(1: trem)]), type="l", ylim=c(0, K*1.3), bty="l",
        ylab="Population size N", xlab='')
    abline(h= K, col="grey50", lty=5, lwd=2)
    abline(v= Tini- trem, col="lightgreen", lty=5, lwd=2)
    lines(smooth.spline(v.N[-(1: trem)], df=10), col="#063d79", lwd=2)
    
    plot(na.omit(v.opt[-(1: trem)]), ylim=c(-1, 8.5),
        type="l", bty='l', xlab='', ylab=expression(paste('Mean optimum ', bar(h))), xlim = c(0,400))
    abline(v= Tini- trem, col="lightgreen", lty=5, lwd=2)
    lines(smooth.spline(na.omit(v.opt[-(1: trem)]), df=10), col="#063d79", lwd=2)  
     
    plot(na.omit(v.tol[-(1:trem)]), ylim=c(0, 8.5),
                 type="l", bty='l', xlab='', ylab=expression(paste('Mean tolerance ', bar(g))), xlim = c(0,400))
    abline(v=Tini-trem, col="lightgreen", lty=5, lwd=2)
    lines(smooth.spline(na.omit(v.tol[-(1:trem)]), df=10), col="#063d79", lwd=2)
    mtext('Generation', 1, outer=T, line=1)
  }) #ends render plot
  
  }) #ends observe event
} #ends server

shinyApp(ui=ui, server=server) #actually creates app!
