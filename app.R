
#install.packages('rJava', type='binary')
library(shiny)
library(ggplot2)
library(shinydashboard)
library(shinythemes)
library(d3heatmap)
#options(java.parameters = "-Xmx5g")
library(rJava)
library(bartMachine)
source("func.R")
source("introduction.R")
source("fun_survbart.R")


actionButton <- function(inputId, label,icon, style = "" ,additionalClass = "") {
    
    if (style %in% c("primary","info","success","warning","danger","inverse","link")) {
        class.style <- paste("btn",style,sep="-")
    } else class.style = ""
    
    
    tags$button(id=inputId, type="button", class=paste("btn action-button ",class.style,additionalClass), label)
}





linearUI <- function(id, label = "linear") {
    ns <- NS(id)
    
    
    header_linear <- dashboardHeader(
        
        title = "Intergrative Bayesian analysis of multi-platform genomics data",
        titleWidth = 700
    )
    
    sidebar_linear <- dashboardSidebar(
        width = 200,
        sidebarMenu(
            id=ns("tabs"),
            menuItem("Introduction", tabName = ns("iintro"), icon = icon("upload", lib = "glyphicon")),
            menuItem("Data Input", tabName = ns("idata"), icon = icon("upload", lib = "glyphicon")),
            menuItem("Posterior Probability Plots", tabName = ns("rdata"), icon = icon("stats", lib = "glyphicon")),
            menuItem("Gene Tables", tabName = ns("gdata"), icon = icon("th", lib = "glyphicon")),
            menuItem("Mechanistic Model Fits", tabName = ns("mcodata"), icon = icon("th", lib = "glyphicon"))
            
        )
    )
    
    body_linear <- dashboardBody(
        tags$head(
            tags$style(HTML("
                    .shiny-progress-container{
                    }

                    "))
        ),
        
        tabItems(
            
            tabItem(tabName = ns("iintro"),
                    
                    fluidRow(img(src='ibag.png',align="left"),
                             
                             p("iBAG is short for integrative Bayesian analysis of high-dimensional multiplatform genomics data. 
                               iBAG is a general framework for integrating information across genomic, transcriptomic and epigenetic 
                               data. Briefly, iBAG uses a novel hierarchical procedure by breaking the modeling into two parts, a 
                               mechanistic component that clarifies the molecular behaviors, mechanisms and relationships between 
                               and within the different types of molecular platforms. Subsequently, a clinical component that utilizes 
                               this information to assess associations between the phenotypes and clinical outcomes that characterize 
                               cancer development and progression (e.g. survival times, treatment arms, response to chemotherapy and tumor 
                               [sub]-types). The  Figure shows a schematic representation of the iBAG modeling strategy. The statistical 
                               formulation of the iBAG models can be found",a('here', href='http://www.ncbi.nlm.nih.gov/pubmed/23142963'),
                               "and", a('here', href='http://www.ncbi.nlm.nih.gov/pubmed/24053265'),". 
                               A standalone version of this code along with an example dataset is available at" , 
                               a('here', href='http://odin.mdacc.tmc.edu/~vbaladan/Veera_Home_Page/iBAG.zip'),"."),
                             
                             box(title = "Input Data", status = "primary",width = 4, background = "black", 
                                 "If you understand the data format to input for the app click the button to proceed", 
                                 actionButton(ns("inButton"), "Input data for iBAG",icon = icon("play"),style="success"))
                             
                    ),
                    fluidRow(
                        
                        box(title = "This app requires the following data in the format 
                            described below", status = "primary",width = 12, background = "black",
                            
                            fluidRow( box(title = "Input format for clinical data", status = "primary",width = 4, background = "black",
                                          p("The rows should represent the individuals and the  2 columns should be the 
                                            individual Id and the response (Binary, Continuous, Survival)."), img(src='sdata.png') ),
                                      
                                      fluidRow(   box(title = "Input format for mRNA data", status = "primary",width = 7, background = "black", 
                                                      p("The rows should be individuals and the columns should be the genes."), img(src='mrnadata.png') ),
                                                  
                                                  
                                                  
                                                  box(title = "Input format for genomic data(Methylation/CNA)", status = "primary",width = 7, background = "black", 
                                                      p("The rows should be individuals and the columns should be markers. Each marker 
                                                        name should include the name of the gene it is present in."), img(src='gdata.png') )
                                                  
                                                  
                                      )
                                      
                            )
                            
                        )
                        
                        
                    )
                    
                    
                    
                    
            ),
            
            
            
            
            tabItem(tabName = ns("idata"),
                    fluidRow(
                        box(title = "Input multi-platform Genomics Data", status = "primary", solidHeader = TRUE, width=6,
                            checkboxInput(ns("cn"), label = "Copy Number"),
                            conditionalPanel(
                                condition = paste0("input.",ns("cn"),"==true"), # "input.cn==true"
                                fileInput(ns('cnfile'), 'Input Copy Number Data File')
                            ),
                            
                            checkboxInput(ns("meth"), label = "DNA Methylation"),
                            conditionalPanel(
                                condition = paste0("input.",ns("meth"),"==true"), # "input.meth==true"
                                fileInput(ns('methfile'), 'Input DNA Methylation Data File')
                            ),
                            checkboxInput(ns("mrna"), label = "MRNA Expression"),
                            conditionalPanel(
                                condition = paste0("input.",ns("mrna"),"==true"), # "input.mrna==true"
                                fileInput(ns('mrnafile'), 'Input mRNA Expression Data File')
                            )),
                        
                        box(title = "Input Clinical Response Data", status = "primary", solidHeader = TRUE, width=6,
                            selectInput(ns("rdata"), label = h3(" Select data type"),
                                        choices = list("None"=0,"Continuous" = 2,
                                                       "Survival (Uncensored)" = 3),selected=0),
                            conditionalPanel(
                                condition = paste0("input.",ns("rdata"),">0"), # "input.rdata>0"
                                fileInput(ns('rfile'), paste("Input clinical data file"))
                            ))),
                    
                    sliderInput(ns("mruns"), "Number of MCMC Runs (Burn in is 5% of runs) 
                                (Calibrate the number of runs according to the computer it is being run on) :",
                                min = 10, max = 10000, value = 10, step= 100),
                    
                    fluidRow(
                        box(title = "Input Data summary", status = "primary",width = 4, background = "black",
                            tableOutput(ns('mytable'))),
                        
                        
                        box(title = "Run Analysis", status = "primary",width = 4, background = "black", 
                            "If you conform with the input data summary press the button to run the analysis", 
                            actionButton(ns("goButton"), "Run iBAG!",icon = icon("play"),style="success"))
                    )
                    
                    
                    
                    
            ),
            
            tabItem(tabName = ns("rdata"),
                    sliderInput(ns("delta"), "Delta:",
                                min = 0, max = 1, value = 0.05, step= 0.01),
                    plotOutput(ns("plot1")),
                    plotOutput(ns("plot2"))
            ),
            tabItem(tabName = ns("gdata"),
                    fluidRow(
                        box(title = "Gene Summary for positive effect on outcome", background = "black",
                            tableOutput(ns("table1"))),
                        box(title = "Gene Summary for negative effect on outcome", background = "black",
                            tableOutput(ns("table2")))
                    )
            ),
            
            tabItem(tabName = ns("mcodata"),
                    fluidRow(
                        uiOutput(ns("col")),
                        box(title = "Mechanistic Model Gene Summary",width=12, background = "black",
                            tableOutput(ns("table3"))),
                        selectInput(ns("palette"), "Choose Color Theme", c("YlOrRd", "RdYlBu", "Greens", "Blues")),
                        d3heatmapOutput(ns("heatmap"),      width="90%",  height="1000px")
                    )
            )
        )
        
        
        
        
        
    )
    
    
    
    dashboardPage(skin = "blue", header_linear, sidebar_linear, body_linear)
}

linearServer <- function(input, output, session) {
    df<-eventReactive(input$goButton, {

        withProgress(message = 'Fitting Mechanistic Model', value = 1/10, {
            nruns<-input$mruns
            GBM_data <- mechmodel(methdata()$data,mrnadata()$data,cndata()$data,rdata()$data)
            incProgress(1/10, message = paste("Fitting Clinical Model"))
            to_gibbs <- prep_and_get_dims(X=GBM_data$X, clinical_response=GBM_data$OurSurvival, take_log=TRUE, GBM=TRUE)
            incProgress(1/10, detail = paste("Initializing MCMC"))
            initial <- get_starting_values_NG(S=nruns, p=to_gibbs$p, k=to_gibbs$k, n=to_gibbs$n, X=to_gibbs$X, Y=to_gibbs$Y,names_to_keep = to_gibbs$names_to_keep)
            incProgress(1/10, detail = paste("Running MCMC"))
            M <- mean(coef(lm(to_gibbs$Y~to_gibbs$X - 1))^2)	#b_tilde=M per Griffin & Brown (2009)
            final <- MC_samples_NG_no_sig_sq_in_beta_prior(PARAM=initial$PARAM, X=to_gibbs$X, Y=to_gibbs$Y, p=to_gibbs$p, 
                                                           k=to_gibbs$k, n=to_gibbs$n,a=0.001, b=0.001, c=1, a_tilde=2, b_tilde=M, 
                                                           tune=0.6, beta_names=initial$beta_names,gam_n2_names=initial$gam_n2_names, 
                                                           lam_names=initial$lam_names, psi_names=initial$psi_names)

            burn_in <- floor(0.05*nruns)+1
            post_means <- apply(final$PARAM[(burn_in+1):nrow(final$PARAM),],2,mean)

            n<-dim(GBM_data$OurMRNA)[1]
            p<-dim(GBM_data$OurMRNA)[2]

            oo <- matrix(nrow=3,ncol=p)
            for (i in 1:(p-1)){ oo[,i] <- c(i,i+p,i+2*p+1) }
            incProgress(6/10, message = paste("iBag Analysis Complete"))
            status<-"iBAG Analysis Complete"
        })
        return(list(to_gibbs=to_gibbs,GBM_data=GBM_data,oo=oo,final=final,initial=initial,status=status,burn_in=burn_in))


    })




    observeEvent(input$goButton, {
        # newtab <- switch(input$tabs,
        #                  "idata" = "rdata",
        #                  "rdata" = "idata"
        # )
        #updateTabItems(session = parent_session, ns("tabs"), newtab)
    }

    )

    observeEvent(input$inButton, {
        newtab <- switch(input$tabs,
                         "iintro" = "idata",
                         "idata" = "iintro"
        )
        updateTabItems(session, "tabs", newtab)
    }

    )


    plotdata<-reactive({
        initial<-df()$initial
        final<-df()$final
        oo<-df()$oo
        burn_in=df()$burn_in
        to_gibbs=df()$to_gibbs
        GBM_data=df()$GBM_data


        delta<-input$delta
        delta_star <- c( log(1-delta), log(1+delta) )
        pos <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t>delta_star[2]))
        neg <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t<delta_star[1]))
        pgene<-which(pos>0.5)
        ngene<-which(neg>0.5)
        cname<- colnames(GBM_data$OurMRNA)
        p<-length(cname)
        pgene1<-pgene[pgene<=p]
        ngene1<-ngene[ngene<=p]
        pgene2<-pgene[(pgene<=2*p) & (pgene>=p+1)]-p
        ngene2<-ngene[(ngene<=2*p) & (ngene>=p+1)]-p
        pgene3<-pgene[(pgene<=3*p) & (pgene>=2*p+1)]-2*p
        ngene3<-ngene[(ngene<=3*p) & (ngene>=2*p+1)]-2*p
        tmp<-array(0,dim=p)
        tmp1<- tmp
        tmp1[pgene1]<-1
        tmp2<- tmp
        tmp2[ pgene2 ]<-1
        tmp3<- tmp
        tmp3[ pgene3]<-1


        t1<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
        maxl<- max(colSums(t1))
        v1<-array("",dim=c(maxl,3))
        c1<-cname[which(t1$meth>0)]
        if(length(c1)>0) v1[1:length(c1),1]<-c1
        c2<-cname[which(t1$cnv>0)]
        if(length(c2)>0) v1[1:length(c2),2]<-c2
        c3<-cname[which(t1$other>0)]
        if(length(c3)>0) v1[1:length(c3),3]<-c3
        colnames(v1)<-c("Methylation","Copy Number", "Other")

        tmp1<- tmp
        tmp1[ngene1]<-1
        tmp2<- tmp
        tmp2[ ngene2 ]<-1
        tmp3<- tmp
        tmp3[ ngene3]<-1


        t2<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
        maxl<- max(colSums(t2))
        v2<-array("",dim=c(maxl,3))
        c1<-cname[which(t2$meth>0)]
        if(length(c1)>0) v2[1:length(c1),1]<-c1
        c2<-cname[which(t2$cnv>0)]
        if(length(c2)>0) v2[1:length(c2),2]<-c2
        c3<-cname[which(t2$other>0)]
        if(length(c3)>0) v2[1:length(c3),3]<-c3
        colnames(v2)<-c("Methylation","Copy Number", "Other")

        aa<-data.frame(x=1:length(initial$beta_names),y=pos[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
        p1<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+
            scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+
            theme(legend.title=element_blank())+ylab("Pr(beta > log(1+delta))")+ggtitle("Posterior Probabilities (Positive)")+xlab("Genes")



        aa<-data.frame(x=1:length(initial$beta_names),y=neg[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
        p2<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+
            scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+
            theme(legend.title=element_blank())+ylab("Pr(beta > log(1-delta))")+ggtitle("Posterior Probabilities (Negative)")+xlab("Genes")

        return(list(p1=p1,p2=p2,t1=v1,t2=v2))
    })

    output$col <- renderUI({

        selectInput("gene", "Select the Gene",  colnames(df()$GBM_data$OurMRNA))
    })

    output$plot1 <- renderPlot({

        plotdata()$p1


    })

    output$plot2 <- renderPlot({

        plotdata()$p2


    })


    output$table1 = renderTable({
        plotdata()$t1
    },
    include.rownames=FALSE)
    output$table2 = renderTable({
        plotdata()$t2
    },
    include.rownames=FALSE)


    output$table3 = renderTable({
        GBM_data=df()$GBM_data
        cname<- colnames(GBM_data$OurMRNA)
        i<-1
        i<- which(cname==input$gene)
        #     tmp<-data.frame(Gene= input$gene,Total =GBM_data$SST[i],Methylation =GBM_data$SSM[i],CopyNumber =GBM_data$SSCN[i] ,Other =GBM_data$SSE[i] )
        #     tmp1<-data.frame(Gene= input$gene,Total =GBM_data$SST[i]/GBM_data$SST[i],Methylation =GBM_data$SSM[i]/GBM_data$SST[i],CopyNumber =GBM_data$SSCN[i]/GBM_data$SST[i] ,Other =GBM_data$SSE[i]/GBM_data$SST[i] )
        #     tmp<-rbind(tmp,tmp1)
        tmp<-c(cname[i], paste(round(GBM_data$SSM[i]/GBM_data$SST[i]*100,2),"%"),paste(round(GBM_data$SSCN[i]/GBM_data$SST[i] *100,2),"%"),paste(round(GBM_data$SSE[i]/GBM_data$SST[i]*100,2),"%") )
        tmp<-data.frame(Gene= tmp[1],Methylation = tmp[2],CopyNumber =tmp[3],Other =tmp[4] )
        return(tmp)
    },
    include.rownames=FALSE)



    output$heatmap <- renderD3heatmap({
        GBM_data=df()$GBM_data
        cname<- colnames(GBM_data$OurMRNA)
        tmp<-data.frame(Methylation =GBM_data$SSM/GBM_data$SST*100,CopyNumber =GBM_data$SSCN/GBM_data$SST *100,Other =GBM_data$SSE/GBM_data$SST*100 ,row.names = cname )
        d3heatmap(
            tmp,
            # scale(mtcars),
            colors = input$palette,
            dendrogram =  "row",
            xaxis_height = 200, yaxis_width = 80,
            xaxis_font_size = 25, yaxis_font_size = 20
        )
    })


    methdata<-reactive({
        methfile <- input$methfile
        if (is.null(methfile)){
            return(NULL)
        }else{
            mdata<-read.csv(methfile$datapath)
            tb=data.frame("Data Type"="Methylation", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }

    })

    cndata<-reactive({
        cnfile <- input$cnfile
        if (is.null(cnfile)){
            return(NULL)
        }else{
            mdata<-read.csv(cnfile$datapath)
            tb=data.frame("Data Type"="Copy Number", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }

    })

    mrnadata<-reactive({
        mrnafile <- input$mrnafile
        if (is.null(mrnafile)){
            return(NULL)
        }else{
            mdata<-read.csv(mrnafile$datapath)
            tb=data.frame("Data Type"="mRNA", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }

    })

    rdata<-reactive({
        rfile <- input$rfile
        if (is.null(rfile)){
            return(NULL)
        }else{
            mdata<-read.csv(rfile$datapath)
            tb=data.frame("Data Type"="Response", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }

    })

    output$mytable = renderTable({
        rbind(methdata()$tb,cndata()$tb,mrnadata()$tb,rdata()$tb)
    },
    include.rownames=FALSE
    )
}


nonlinearUI <- function(id, label = "nonlinear") {
    ns <- NS(id)
    
    
    header_nonlinear <- dashboardHeader(
        
        title = "Intergrative Bayesian analysis of multi-platform genomics data",
        titleWidth = 700
    )
    
    sidebar_nonlinear <- dashboardSidebar(
        width = 200,
        sidebarMenu(
            id=ns("tabs"),
            menuItem("Introduction", tabName = ns("iintro"), icon = icon("upload", lib = "glyphicon")),
            menuItem("Data Input", tabName = ns("idata"), icon = icon("upload", lib = "glyphicon")),
            menuItem("Mechanistic Model Fits", tabName = ns("mcodata"), icon = icon("stats", lib = "glyphicon")),
            menuItem("Posterior Probability Plots", tabName = ns("rdata"), icon = icon("stats", lib = "glyphicon")),
            menuItem("Gene Tables", tabName = ns("gdata"), icon = icon("th", lib = "glyphicon"))
        )
    )
    
    
    
    body_nonlinear <- dashboardBody(
        tags$head(
            tags$style(HTML("
                    .shiny-progress-container{
                    }
                    
                    "))
        ),
        
        tabItems(
            
            tabItem(tabName = ns("iintro"),
                    
                    fluidRow(img(src='ibag.png',align="left"),
                             
                             p("iBAG is short for integrative Bayesian analysis of high-dimensional multiplatform genomics data. iBAG is a general 
                               framework for integrating information across genomic, transcriptomic and epigenetic data. Briefly, iBAG uses a novel 
                               hierarchical procedure by breaking the modeling into two parts, a mechanistic component that clarifies the molecular 
                               behaviors, mechanisms and relationships between and within the different types of molecular platforms. Subsequently, 
                               a clinical component that utilizes this information to assess associations between the phenotypes and clinical outcomes 
                               that characterize cancer development and progression (e.g. survival times, treatment arms, response to chemotherapy and 
                               tumor [sub]-types). The  Figure shows a schematic representation of the iBAG modeling strategy. The statistical formulation 
                               of the iBAG models can be found",a('here', href='http://www.ncbi.nlm.nih.gov/pubmed/23142963'),"and", 
                               a('here', href='http://www.ncbi.nlm.nih.gov/pubmed/24053265'),". A standalone version of this code along with an example dataset is available at" , 
                               a('here', href='http://odin.mdacc.tmc.edu/~vbaladan/Veera_Home_Page/iBAG.zip'),"."),
                             box(title = "Input Data", status = "primary",width = 4, background = "black", 
                                 "If you understand the data format to input for the app click the button to proceed", 
                                 actionButton(ns("inButton"), "Input data for iBAG",icon = icon("play"),style="success"))
                             
                    ),
                    fluidRow(
                        
                        
                        
                        box(title = "This app requires the following data in the format described below", 
                            status = "primary",width = 12, background = "black",
                            
                            fluidRow( box(title = "Input format for clinical data", status = "primary",width = 4, 
                                          background = "black",p("The rows should represent the individuals and the  2 columns should be the individual 
                                                                 Id and the response (Binary, Continuous, Survival)."), img(src='sdata.png') ),
                                      
                                      fluidRow(   box(title = "Input format for mRNA data", status = "primary",width = 7, 
                                                      background = "black", p("The rows should be individuals and the columns should be the genes."), img(src='mrnadata.png') ),
                                                  
                                                  
                                                  
                                                  box(title = "Input format for genomic data(Methylation/CNA)", status = "primary",width = 7, background = "black", 
                                                      p("The rows should be individuals and the columns should be markers. Each marker name should include the name 
                                                        of the gene it is present in."), img(src='gdata.png') )
                                                  
                                                  
                                      )
                                      
                            )
                            
                        )
                        
                        
                    )
                    
                    
                    
                    
            ),
            
            
            
            
            tabItem(tabName = ns("idata"),
                    fluidRow(
                        box(title = "Input multi-platform Genomics Data", status = "primary", solidHeader = TRUE, width=6,
                            checkboxInput(ns("cn"), label = "Copy Number"),
                            conditionalPanel(
                                condition = paste0("input.",ns("cn"),"==true"), # "input.cn==true"
                                fileInput(ns('cnfile'), 'Input Copy Number Data File')
                            ),
                            
                            checkboxInput(ns("meth"), label = "DNA Methylation"),
                            conditionalPanel(
                                condition = paste0("input.",ns("meth"),"==true"), # "input.meth==true"
                                fileInput(ns('methfile'), 'Input DNA Methylation Data File')
                            ),
                            checkboxInput(ns("mrna"), label = "MRNA Expression"),
                            conditionalPanel(
                                condition = paste0("input.",ns("mrna"),"==true"), # "input.mrna==true"
                                fileInput(ns('mrnafile'), 'Input mRNA Expression Data File')
                            )),
                        
                        box(title = "Input Clinical Response Data", status = "primary", solidHeader = TRUE, width=6,
                            selectInput(ns("rdata"), label = h3(" Select data type"),
                                        choices = list("None"=0,"Continuous" = 2,
                                                       "Survival (Uncensored)" = 3),selected=0),
                            conditionalPanel(
                                condition = paste0("input.",ns("rdata"),">0"), # "input.rdata>0"
                                fileInput(ns('rfile'), paste("Input clinical data file"))
                            ))),
                    
                    sliderInput(ns("mruns"), "Number of MCMC Runs (Burn in is 5% of runs) (Calibrate the number of runs according to the computer it is being run on) :",
                                min = 10, max = 10000, value = 10, step= 100),
                    
                    fluidRow(
                        box(title = "Input Data summary", status = "primary",width = 4, background = "black",
                            tableOutput(ns('mytable'))),
                        
                        
                        box(title = "Run Analysis", status = "primary",width = 4, background = "black", "If you conform with the input data summary press the button to run the analysis", 
                            actionButton(ns("goButton1"), "Run linear iBAG!",icon = icon("play"),style="success"), 
                            actionButton(ns("goButton2"), "Run bart iBAG!",icon = icon("play"),style="success")),
                        box(title="R square:", status = "primary",width = 4, background = "black",textOutput("rsq"))
                    )
                    
                    
                    
                    
            ),
            
            tabItem(tabName = ns("rdata"),
                    sliderInput(ns("delta"), "Delta:",
                                min = 0, max = 0.5, value = 0.01, step= 0.001),
                    plotOutput(ns("plot1")),
                    plotOutput(ns("plot2"))
            ),
            tabItem(tabName = ns("gdata"),
                    fluidRow(
                        box(title = "Gene Summary for positive effect on outcome", background = "black",
                            tableOutput(ns("table1"))),
                        box(title = "Gene Summary for negative effect on outcome", background = "black",
                            tableOutput(ns("table2")))
                    )
            ),
            
            tabItem(tabName = ns("mcodata"),
                    fluidRow(
                        uiOutput(ns("col")),
                        box(title = "Mechanistic Model Gene Summary",width=12, background = "black",
                            actionButton(ns("goButton3"), "Go!",icon = icon("play"),style="success"),
                            tableOutput(ns("table3"))),
                        selectInput(ns("palette"), "Choose Color Theme", c("YlOrRd", "RdYlBu", "Greens", "Blues")),
                        d3heatmapOutput(ns("heatmap"),      width="90%",  height="1000px")
                    )
            )
        )
        
        
        
        
        
    )
    
    
    
    
    
    dashboardPage(skin = "blue",header_nonlinear,sidebar_nonlinear,body_nonlinear)
}

nonlinearServer <- function(input, output, session) {
    
    global <- reactiveValues()
    
    df<-eventReactive(input$goButton2, {
        
        withProgress(message = 'Fitting Mechanistic Model', value = 1/10, {
            nruns<-input$mruns
            GBM_data <- mechmodel(methdata()$data,mrnadata()$data,cndata()$data,rdata()$data)
            #GBM_data <-mechmodel(meth,mrna,cnv,dsurv)
            incProgress(1/10, message = paste("Fitting Clinical Model"))
            
            to_gibbs <- prep_and_get_dims(X=GBM_data$X, clinical_response = GBM_data$OurSurvival,take_log=TRUE, GBM=TRUE)
            
            #initial <- get_starting_values_NG(S=nruns, p=to_gibbs$p, k=to_gibbs$k, n=to_gibbs$n, X=to_gibbs$X, Y=to_gibbs$Y,names_to_keep = to_gibbs$names_to_keep)
            #incProgress(1/10, detail = paste("Running MCMC"))
            #M <- mean(coef(lm(to_gibbs$Y~to_gibbs$X - 1))^2)	#b_tilde=M per Griffin & Brown (2009)
            #final <- MC_samples_NG_no_sig_sq_in_beta_prior(PARAM=initial$PARAM, X=to_gibbs$X, Y=to_gibbs$Y, p=to_gibbs$p, k=to_gibbs$k, n=to_gibbs$n,a=0.001, b=0.001, c=1, a_tilde=2, b_tilde=M, tune=0.6, beta_names=initial$beta_names,gam_n2_names=initial$gam_n2_names, lam_names=initial$lam_names, psi_names=initial$psi_names)
            
            burn_in <- floor(0.05*nruns)+1
            #post_means <- apply(final$PARAM[(burn_in+1):nrow(final$PARAM),],2,mean)
            
            n<-dim(GBM_data$OurMRNA)[1]
            p<-dim(GBM_data$OurMRNA)[2]
            
            oo <- matrix(nrow=3,ncol=p)
            for (i in 1:(p-1)){ oo[,i] <- c(i,i+p,i+2*p+1) }
            m_c= bartMachine(X=data.frame(to_gibbs$X), y=to_gibbs$Y,num_trees = 100)
            var_imp=investigate_var_importance(m_c, num_replicates_for_avg = 20)
            rsq=m_c$PseudoRsq
            delta<-input$delta
            delta_star <- c( log(1-delta), log(1+delta) )
            pos=var_imp$avg_var_props>delta_star[2]
            neg=var_imp$avg_var_props<delta_star[1]
            
            
            incProgress(6/10, message = paste("bart iBag Analysis for survival data Complete"))
            status<-"bart iBAG Analysis Complete"
            
            
            # pos <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t>delta_star[2]))
            # neg <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t<delta_star[1]))
            aname=colnames(GBM_data$X)
            cname<- colnames(GBM_data$OurMRNA)
            po<-aname%in%names(which(pos))
            pgene<-which(po)
            ne<-aname%in%names(which(neg))
            ngene<-which(ne)
            
            p<-length(cname)
            pgene1<-pgene[pgene<=p]
            ngene1<-ngene[ngene<=p]
            pgene2<-pgene[(pgene<=2*p) & (pgene>=p+1)]-p
            ngene2<-ngene[(ngene<=2*p) & (ngene>=p+1)]-p
            pgene3<-pgene[(pgene<=3*p) & (pgene>=2*p+1)]-2*p
            ngene3<-ngene[(ngene<=3*p) & (ngene>=2*p+1)]-2*p
            tmp<-array(0,dim=p)
            tmp1<- tmp
            tmp1[pgene1]<-1
            tmp2<- tmp
            tmp2[ pgene2 ]<-1
            tmp3<- tmp
            tmp3[ pgene3]<-1
            
            
            t1<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
            maxl<- max(colSums(t1))
            v1<-array("",dim=c(maxl,3))
            c1<-cname[which(t1$meth>0)]
            if(length(c1)>0) v1[1:length(c1),1]<-c1
            c2<-cname[which(t1$cnv>0)]
            if(length(c2)>0) v1[1:length(c2),2]<-c2
            c3<-cname[which(t1$other>0)]
            if(length(c3)>0) v1[1:length(c3),3]<-c3
            colnames(v1)<-c("Methylation","Copy Number", "Other")
            
            tmp1<- tmp
            tmp1[ngene1]<-1
            tmp2<- tmp
            tmp2[ ngene2 ]<-1
            tmp3<- tmp
            tmp3[ ngene3]<-1
            
            
            t2<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
            maxl<- max(colSums(t2))
            v2<-array("",dim=c(maxl,3))
            c1<-cname[which(t2$meth>0)]
            if(length(c1)>0) v2[1:length(c1),1]<-c1
            c2<-cname[which(t2$cnv>0)]
            if(length(c2)>0) v2[1:length(c2),2]<-c2
            c3<-cname[which(t2$other>0)]
            if(length(c3)>0) v2[1:length(c3),3]<-c3
            colnames(v2)<-c("Methylation","Copy Number", "Other")
            
            #aa<-data.frame(x=1:length(initial$beta_names),y=pos[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
            aa<-data.frame(x=1:length(cname),y=var_imp$avg_var_props[aname][oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
            p1<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=delta_star[2]))+theme(legend.title=element_blank())+ylab("Posterior Inclusion Proportion that > log(1+delta))")+ggtitle("Variable Importance")+xlab("Genes")
            
            
            p2<-NULL
        })
        return(list(to_gibbs=to_gibbs,GBM_data=GBM_data,oo=oo,status=status,burn_in=burn_in,p1=p1,p2=p2,t1=v1,t2=v2,rsq=rsq))
        
        
    })
    
    df1<-eventReactive(input$goButton1, {
        
        withProgress(message = 'Fitting Mechanistic Model', value = 1/10, {
            nruns<-input$mruns
            GBM_data <- mechmodel(methdata()$data,mrnadata()$data,cndata()$data,rdata()$data)
            incProgress(1/10, message = paste("Fitting Clinical Model"))
            to_gibbs <- prep_and_get_dims(X=GBM_data$X,clinical_response = GBM_data$OurSurvival,  take_log=TRUE, GBM=TRUE)
            incProgress(1/10, detail = paste("Initializing MCMC"))
            initial <- get_starting_values_NG(S=nruns, p=to_gibbs$p, k=to_gibbs$k, n=to_gibbs$n, X=to_gibbs$X, Y=to_gibbs$Y,names_to_keep = to_gibbs$names_to_keep)
            incProgress(1/10, detail = paste("Running MCMC"))
            M <- mean(coef(lm(to_gibbs$Y~to_gibbs$X - 1))^2)	#b_tilde=M per Griffin & Brown (2009)
            final <- MC_samples_NG_no_sig_sq_in_beta_prior(PARAM=initial$PARAM, X=to_gibbs$X, Y=to_gibbs$Y, p=to_gibbs$p, k=to_gibbs$k, n=to_gibbs$n,a=0.001, b=0.001, c=1, a_tilde=2, b_tilde=M, tune=0.6, beta_names=initial$beta_names,gam_n2_names=initial$gam_n2_names, lam_names=initial$lam_names, psi_names=initial$psi_names)
            
            burn_in <- floor(0.05*nruns)+1
            post_means <- apply(final$PARAM[(burn_in+1):nrow(final$PARAM),],2,mean)
            
            n<-dim(GBM_data$OurMRNA)[1]
            p<-dim(GBM_data$OurMRNA)[2]
            
            oo <- matrix(nrow=3,ncol=p)
            for (i in 1:(p-1)){ oo[,i] <- c(i,i+p,i+2*p+1) }
            incProgress(6/10, message = paste("linear iBag Analysis Complete"))
            status<-"linear iBAG Analysis Complete"
            delta<-input$delta
            delta_star <- c( log(1-delta), log(1+delta) )
            pos <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t>delta_star[2]))
            neg <- apply(final$PARAM[-(1:burn_in),initial$beta_names],2,function(t) mean(t<delta_star[1]))
            pgene<-which(pos>0.5)
            ngene<-which(neg>0.5)
            cname<- colnames(GBM_data$OurMRNA)
            p<-length(cname)
            pgene1<-pgene[pgene<=p]
            ngene1<-ngene[ngene<=p]
            pgene2<-pgene[(pgene<=2*p) & (pgene>=p+1)]-p
            ngene2<-ngene[(ngene<=2*p) & (ngene>=p+1)]-p
            pgene3<-pgene[(pgene<=3*p) & (pgene>=2*p+1)]-2*p
            ngene3<-ngene[(ngene<=3*p) & (ngene>=2*p+1)]-2*p
            tmp<-array(0,dim=p)
            tmp1<- tmp
            tmp1[pgene1]<-1
            tmp2<- tmp
            tmp2[ pgene2 ]<-1
            tmp3<- tmp
            tmp3[ pgene3]<-1
            
            
            t1<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
            maxl<- max(colSums(t1))
            v1<-array("",dim=c(maxl,3))
            c1<-cname[which(t1$meth>0)]
            if(length(c1)>0) v1[1:length(c1),1]<-c1
            c2<-cname[which(t1$cnv>0)]
            if(length(c2)>0) v1[1:length(c2),2]<-c2
            c3<-cname[which(t1$other>0)]
            if(length(c3)>0) v1[1:length(c3),3]<-c3
            colnames(v1)<-c("Methylation","Copy Number", "Other")
            
            tmp1<- tmp
            tmp1[ngene1]<-1
            tmp2<- tmp
            tmp2[ ngene2 ]<-1
            tmp3<- tmp
            tmp3[ ngene3]<-1
            
            
            t2<-data.frame(meth=tmp1,cnv=tmp2,other=tmp3,row.names = cname)
            maxl<- max(colSums(t2))
            v2<-array("",dim=c(maxl,3))
            c1<-cname[which(t2$meth>0)]
            if(length(c1)>0) v2[1:length(c1),1]<-c1
            c2<-cname[which(t2$cnv>0)]
            if(length(c2)>0) v2[1:length(c2),2]<-c2
            c3<-cname[which(t2$other>0)]
            if(length(c3)>0) v2[1:length(c3),3]<-c3
            colnames(v2)<-c("Methylation","Copy Number", "Other")
            
            aa<-data.frame(x=1:length(initial$beta_names),y=pos[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
            p1<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+theme(legend.title=element_blank())+ylab("Pr(beta > log(1+delta))")+ggtitle("Posterior Probabilities (Positive)")+xlab("Genes")
            
            
            
            aa<-data.frame(x=1:length(initial$beta_names),y=neg[oo],g=c(rep(1,to_gibbs$p[1]),rep(2,to_gibbs$p[2]),rep(3,to_gibbs$p[3]))[oo])
            p2<-ggplot()+geom_bar(aes(x=x,y=y,fill = as.factor(g)), position = "dodge", stat="identity",data=aa)+scale_fill_discrete(labels=c("Methylation","Copy Number","Other"))+scale_x_continuous(breaks=c(seq(1,144,by=3)), labels=c(colnames(GBM_data$OurMRNA))  )+theme(axis.text.x = element_text(angle=90))+geom_hline(aes(yintercept=0.5))+theme(legend.title=element_blank())+ylab("Pr(beta > log(1-delta))")+ggtitle("Posterior Probabilities (Negative)")+xlab("Genes")
            
        })
        return(list(to_gibbs=to_gibbs,GBM_data=GBM_data,oo=oo,final=final,initial=initial,status=status,burn_in=burn_in,p1=p1,p2=p2,t1=v1,t2=v2))
        
        
    })
    
    
    
    
    observeEvent(input$goButton1, {
        newtab <- switch(input$tabs,
                         "idata" = "mcodata",
                         "mcodata" = "idata"
        )
        updateTabItems(session, ns("tabs"), newtab)
        global$p1=df1()$p1
        global$p2=df1()$p2
        global$GBM_data=df1()$GBM_data
        global$cname=colnames(df1()$GBM_data$OurMRNA)
    }
    
    )
    
    observeEvent(input$goButton2, {
        newtab <- switch(input$tabs,
                         "idata" = "mcodata",
                         "mcodata" = "idata"
        )
        updateTabItems(session, ns("tabs"), newtab)
        global$p1=df()$p1
        global$p2=df()$p2
        global$GBM_data=df()$GBM_data
        global$cname=colnames(df()$GBM_data$OurMRNA)
    }
    
    )
    
    
    observeEvent(input$goButton3, {
        i=1
        cname=global$cname
        i<- which(cname==input$gene)
        GBM_data=global$GBM_data
        tmp<-c(cname[i], paste(round(GBM_data$SSM[i]/GBM_data$SST[i]*100,2),"%"),paste(round(GBM_data$SSCN[i]/GBM_data$SST[i] *100,2),"%"),paste(round(GBM_data$SSE[i]/GBM_data$SST[i]*100,2),"%") )
        global$table3<-data.frame(Gene= tmp[1],Methylation = tmp[2],CopyNumber =tmp[3],Other =tmp[4] )
    })
    
    
    observeEvent(input$inButton, {
        newtab <- switch(input$tabs,
                         "iintro" = "idata",
                         "idata" = "iintro"
        )
        updateTabItems(session, ns("tabs"), newtab)
    }
    
    )
    
    
    
    output$col <- renderUI({
        
        selectInput(ns("gene"), "Select the Gene",  global$cname)
    })
    
    
    output$rsq<-renderText({
        df()$rsq
    })
    
    output$plot1 <- renderPlot({
        
        global$p1
        
        
    })
    
    
    
    output$plot2 <- renderPlot({
        
        
        global$p2
        
    })
    
    
    
    #  data <- reactive({
    #    validate(
    #      need(input$gene2 != "", "Please select a gene")
    #    )
    #  })
    
    
    
    
    output$table1 = renderTable({
        df()$t1
    },
    include.rownames=FALSE)
    
    output$table2 = renderTable({
        df()$t2
    },
    include.rownames=FALSE)
    
    
    output$table3 = renderTable({
        global$table3
    },
    include.rownames=FALSE)
    
    
    
    
    output$heatmap <- renderD3heatmap({
        GBM_data=global$GBM_data
        cname<- colnames(GBM_data$OurMRNA)
        tmp<-data.frame(Methylation =GBM_data$SSM/GBM_data$SST*100,CopyNumber =GBM_data$SSCN/GBM_data$SST *100,Other =GBM_data$SSE/GBM_data$SST*100 ,row.names = cname )
        d3heatmap(
            tmp,
            # scale(mtcars),
            colors = input$palette,
            dendrogram =  "row",
            xaxis_height = 200, yaxis_width = 80,
            xaxis_font_size = 25, yaxis_font_size = 20
        )
    })
    
    
    methdata<-reactive({
        methfile <- input$methfile
        if (is.null(methfile)){
            return(NULL)
        }else{
            mdata<-read.csv(methfile$datapath)
            tb=data.frame("Data Type"="Methylation", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }
        
    })
    
    cndata<-reactive({
        cnfile <- input$cnfile
        if (is.null(cnfile)){
            return(NULL)
        }else{
            mdata<-read.csv(cnfile$datapath)
            tb=data.frame("Data Type"="Copy Number", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }
        
    })
    
    mrnadata<-reactive({
        mrnafile <- input$mrnafile
        if (is.null(mrnafile)){
            return(NULL)
        }else{
            mdata<-read.csv(mrnafile$datapath)
            tb=data.frame("Data Type"="mRNA", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }
        
    })
    
    rdata<-reactive({
        rfile <- input$rfile
        if (is.null(rfile)){
            return(NULL)
        }else{
            mdata<-read.csv(rfile$datapath)
            tb=data.frame("Data Type"="Response", p=dim(mdata)[2],n=dim(mdata)[1],row.names = NULL)
            return(list(tb=tb,data=mdata))
        }
        
    })
    
    output$mytable = renderTable({
        rbind(methdata()$tb,cndata()$tb,mrnadata()$tb,rdata()$tb)
    },
    include.rownames=FALSE
    )
}


# Define UI for application that draws a histogram
ui <- dashboardPage(
    # Header of the page
    dashboardHeader(

        # # Dashboard title
        # title = "iBAG Dashboard"
        disable = TRUE
    ),

    # Sidebar with all the individual app options
    dashboardSidebar(

        # sidebarMenuOutput(
        #     "MENU_TABS"
        # )
        
        sidebarMenu(
            menuItem("iBag Introduction", tabName = "intro", icon = icon("dashboard")),
            menuItem("Linear iBag", icon = icon("th"), tabName = "input-linear",
                     badgeLabel = "new", badgeColor = "green"),
            menuItem("Non-linear iBag", icon = icon("th"), tabName = "input-nonlinear",
                     badgeLabel = "new", badgeColor = "green")
        )

    ),

    # Body: intro page and each respective app
    dashboardBody(
        
        # A list of all tab items for main page content
        tabItems(
            
            # Introduction tab content to be displayed The contents introduction_page is stored
            # in a separate file so this page does not get cluttered                             # introduction_page
            tabItem(tabName = "intro",
                    introduction_page
            ),
            
            # Linear iBAG tab -------------------------------------------------
            tabItem(tabName = "input-linear",
                linearUI("linear1", "linear 1")
            ),
            # 
            tabItem(tabName = "input-nonlinear",
                    nonlinearUI("nonlinear1", "nonlinear 1")
            )
            
            
        )
        

    )
)


#ui <- linearUI("linear", "linear1")


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    #
    callModule(nonlinearServer, "nonlinear1")
    callModule(linearServer, "linear1")
}

# Run the application 
shinyApp(ui = ui, server = server)
