
DIFshiny<-function(){

ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title = "DIFFERENTIAL ITEM FUNCTIONING",titleWidth = 500),
                    dashboardSidebar(
                      sidebarMenu(id="sidebarmenu",
                                  menuItem("About Package",tabName = "info",icon = icon("info-circle")),
                                  menuItem("DIF Detection",tabName = "diff",icon = icon("cogs")),
                                  menuItem("DIF Simulation",tabName = "dif",icon = icon("cogs")),
                                  menuItem("Data Preview",tabName = "preview",icon = icon("table"))


                      )),
                    dashboardBody(
                      tabItems(

                        tabItem(
                          tabName = "diff",
                          fluidRow(column(8,
                                          box(title="DIF Outputs",verbatimTextOutput("difreel"),solidHeader = TRUE,status = "info",width = 12)),
                                   column(4,
                                          box(title = "UpLoad File",status = "warning",solidHeader = TRUE,width = 12,
                                              fileInput("dataset2", "Please upload your data file in one of .xls, .xslx or .csv formats. Last column of data must be two category group variable. Note: Code 1 for reference group members, 2 for focal group members", placeholder="File",buttonLabel = "Browse",accept = c("xlsx","xls","csv")
                                              )),
                                          box(title = "Method Selection",status = "success",solidHeader = TRUE,width = 12, uiOutput(outputId="methodselect2"),submitButton("Submit")))
                          )),

                        tabItem(tabName = "info",
                                fluidRow(

                                  box(title = "About Package", solidHeader = TRUE, status = "info",htmlOutput("info"),width = 12)
                                )),
                        tabItem(tabName = "preview",
                                fluidRow(

                                  box(title = "Data Preview", solidHeader = TRUE, status = "info",dataTableOutput("preview"),width = 12)
                                )),

                        tabItem(tabName = "dif",
                                fluidRow(column(8,

                                                box(title="DMF Output",verbatimTextOutput("difoutput"),solidHeader = TRUE,status = "info",width = 12)),
                                         column(4,
                                                box(title="Method Selection",solidHeader = TRUE,status = "warning",width = 12,
                                                    uiOutput(outputId="methodselect"),submitButton("Submit")),
                                                box(title = "Select Conditions",status = "warning",solidHeader = TRUE,width = 12,
                                                    uiOutput(outputId="aparmean"),
                                                    uiOutput(outputId="aparsd"),
                                                    uiOutput(outputId="bparmin"),
                                                    uiOutput(outputId="bparmax"),
                                                    uiOutput(outputId="thetamean"),
                                                    uiOutput(outputId="thetasd"),
                                                    uiOutput(outputId="itemnumber"),
                                                    uiOutput(outputId="samplesize"),
                                                    uiOutput(outputId="difratio"),
                                                    uiOutput(outputId="difsize"),submitButton("Submit")

                                                ))),



                        ))

))


server <- function(input, output,session) {


  mydifdata<-reactive({
    inFile <- input$dataset2
    if (is.null(inFile))
      return("Please upload data")
    dataset<- read_xlsx(inFile$datapath, sheet=1)
    data<-as.data.frame(dataset)
    data
  })



  output$methodselect<-renderUI({
    selectInput(inputId ="methodselect",label = NULL,choices = c("Mantel Haenszel","Logistic Regression","Lord","SIBTEST","Raju"))
  })

  output$methodselect2<-renderUI({
    selectInput(inputId ="methodselect2",label = NULL,choices = c("Mantel Haenszel","Logistic Regression","Lord","SIBTEST","Raju"))
  })



  output$aparmean<-renderUI({
    sliderInput(inputId ="aparmean",label = "Mean of a parameters",min = 0.5,max = 1.5, step = 0.1,value = 1)
  })

  output$aparsd<-renderUI({
    sliderInput(inputId ="aparsd",label = "Sd of a parameters",min = 0.05,max = 0.15, step = 0.01,value = 0.01)
  })

  output$bparmin<-renderUI({
    sliderInput(inputId ="bparmin",label = "min value of b parameters",min = -3,max = -1, step = 0.5,value = -2)
  })

  output$bparmax<-renderUI({
    sliderInput(inputId ="bparmax",label = "max value of b parameters",min = 1,max = 3, step = 0.5,value = 2)
  })

  output$thetamean<-renderUI({
    sliderInput(inputId ="thetamean",label = "Mean Theta",min = -0.5,max = 0.5, step = 0.1,value = 0)
  })

  output$thetasd<-renderUI({
    sliderInput(inputId ="thetasd",label = "Sd Theta",min = 0.5,max = 1.5, step = 0.1,value = 1)
  })

  output$itemnumber<-renderUI({
    sliderInput(inputId ="itemnumber",label = "Number of items",min = 20,max = 100, step = 10,value = 20)
  })

  output$samplesize<-renderUI({
    sliderInput(inputId ="samplesize",label = "Sample size",min = 250,max = 2000, step = 250,value = 500)
  })

  output$difsize<-renderUI({
    sliderInput(inputId ="difsize",label = "DIF magnitude",min = 0.25,max = 1.5, step = 0.25,value = 0.50)
  })

  output$difratio<-renderUI({
    sliderInput(inputId ="difratio",label = "Ratio of items have DIF",min = 0.1,max = 0.3, step = 0.1,value = 0.1)
  })

  simdata<-reactive({

    pij <- function(a,b,theta) {1/(1+exp(-a*1.702*(theta-b)))}
    a <- rnorm(input$itemnumber, input$aparmean, input$aparsd)

    b <- runif(input$itemnumber, input$bparmin, input$bparmax)

    refresponse <- matrix(0, input$samplesize, input$itemnumber)



    theta <- rnorm(input$samplesize, input$thetamean,input$thetasd)


    for( i in 1:input$samplesize) {
      for( j in 1:input$itemnumber ) {
        refresponse[i,j]<-ifelse(pij(a=a[j], b=b[j], theta[i]) < runif(1) , 0 ,1)
      }
    }

    ################################

    focresponse <- matrix(0, input$samplesize, input$itemnumber)


    theta <- rnorm(input$samplesize, input$thetamean,input$thetasd)

    bf<-b
    bf[1:input$itemnumber*input$difratio]<-bf[1:input$itemnumber*input$difratio]+input$difsize



    for( i in 1:input$samplesize ) {
      for( j in 1:input$itemnumber ) {
        focresponse[i,j]<-ifelse(pij(a=a[j], b=bf[j], theta[i]) < runif(1) , 0 ,1)
      }
    }

    #################################


    newdata<-rbind(focresponse,refresponse)
    gr<-c(rep(1,input$samplesize),rep(2, input$samplesize))
    Response<-cbind(newdata,gr)

    Response

  })

  output$preview<-renderDataTable({

   if(is.null(input$samplesize))
     return(NULL)

    if(!is.null(simdata())){

      data<-simdata()
    }
    return(data)
  })

  simdata<-reactive({

    pij <- function(a,b,theta) {1/(1+exp(-a*1.702*(theta-b)))}
    a <- rnorm(input$itemnumber, input$aparmean, input$aparsd)

    b <- runif(input$itemnumber, input$bparmin, input$bparmax)

    refresponse <- matrix(0, input$samplesize, input$itemnumber)



    theta <- rnorm(input$samplesize, input$thetamean,input$thetasd)


    for( i in 1:input$samplesize) {
      for( j in 1:input$itemnumber ) {
        refresponse[i,j]<-ifelse(pij(a=a[j], b=b[j], theta[i]) < runif(1) , 0 ,1)
      }
    }

    ################################

    focresponse <- matrix(0, input$samplesize, input$itemnumber)


    theta <- rnorm(input$samplesize, input$thetamean,input$thetasd)

    bf<-b
    bf[1:input$itemnumber*input$difratio]<-bf[1:input$itemnumber*input$difratio]+input$difsize



    for( i in 1:input$samplesize ) {
      for( j in 1:input$itemnumber ) {
        focresponse[i,j]<-ifelse(pij(a=a[j], b=bf[j], theta[i]) < runif(1) , 0 ,1)
      }
    }

    #################################


    newdata<-rbind(focresponse,refresponse)
    gr<-c(rep(1,input$samplesize),rep(2, input$samplesize))
    Response<-cbind(newdata,gr)

    Response

  })



  output$difoutput<-renderPrint({

    if (is.null(input$samplesize))
      return("Please select your simulation conditions and method. ")



    pij <- function(a,b,theta) {1/(1+exp(-a*1.702*(theta-b)))}
    a <- rnorm(input$itemnumber, input$aparmean, input$aparsd)

    b <- runif(input$itemnumber, input$bparmin, input$bparmax)

    refresponse <- matrix(0, input$samplesize, input$itemnumber)



    theta <- rnorm(input$samplesize, input$thetamean,input$thetasd)


    for( i in 1:input$samplesize) {
      for( j in 1:input$itemnumber ) {
        refresponse[i,j]<-ifelse(pij(a=a[j], b=b[j], theta[i]) < runif(1) , 0 ,1)
      }
    }

    ################################

    focresponse <- matrix(0, input$samplesize, input$itemnumber)


    theta <- rnorm(input$samplesize, input$thetamean,input$thetasd)

    bf<-b
    bf[1:input$itemnumber*input$difratio]<-bf[1:input$itemnumber*input$difratio]+input$difsize



    for( i in 1:input$samplesize ) {
      for( j in 1:input$itemnumber ) {
        focresponse[i,j]<-ifelse(pij(a=a[j], b=bf[j], theta[i]) < runif(1) , 0 ,1)
      }
    }

    #################################


    newdata<-rbind(focresponse,refresponse)
    gr<-c(rep(1,input$samplesize),rep(2, input$samplesize))
    Response<-cbind(newdata,gr)



    if(input$methodselect=="Mantel Haenszel"){
      datasayi<-ncol(Response)-1
      sonuc<-difMH(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2)

    }
    if(input$methodselect=="Lord"){
      datasayi<-ncol(Response)-1
      sonuc<-difLord(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2,model="2PL")

    }

    if(input$methodselect=="Raju"){
      datasayi<-ncol(Response)-1
      sonuc<-difRaju(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2,model="2PL")

    }

    if(input$methodselect=="Logistic Regression"){
      datasayi<-ncol(Response)-1
      sonuc<-difLogistic(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2)

    }


    if(input$methodselect=="SIBTEST"){
      datasayi<-ncol(Response)-1
      sonuc<-difSIBTEST(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2)

    }

    sonuc


  })
  output$difreel<-renderPrint({


    Response<- mydifdata()

    if(is.null(input$dataset2))
      return("Please upload your data.")

    if(is.null(input$methodselect2))
      return("Please select your method.")

    if(input$methodselect2=="Mantel Haenszel"){
      datasayi<-ncol(Response)-1
      sonuc<-difMH(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2)

    }
    if(input$methodselect2=="Lord"){
      datasayi<-ncol(Response)-1
      sonuc<-difLord(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2,model="2PL")

    }

    if(input$methodselect2=="Raju"){
      datasayi<-ncol(Response)-1
      sonuc<-difRaju(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2,model="2PL")

    }

    if(input$methodselect2=="Logistic Regression"){
      datasayi<-ncol(Response)-1
      sonuc<-difLogistic(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2)

    }


    if(input$methodselect2=="SIBTEST"){
      datasayi<-ncol(Response)-1
      sonuc<-difSIBTEST(Data=Response[,1:datasayi],group=Response[,ncol(Response)],focal.name=2)

    }

    sonuc


  })

  output$info<-renderText({
    paste(p(strong('Package:'), "DIFshiny"),p(strong('Background Packages:'), "difR,","shiny,","shinydashboard"),
          p(strong('Package Description:'), "This package helps to Differential Item Functioning (DIF) Analysis with shiny application interfaces.    You can run the functions in this package without any arguments and perform your DIF analysis using user-friendly interfaces."),

          p(strong('Package Author:'), "Huseyin Yildiz"),
          p(strong('e-mail:'), tags$a(href="mailto:huseyinyildiz35@gmail.com", "huseyinyildiz35@gmail.com"))

          )
  })



}

shinyApp(ui, server)


}



