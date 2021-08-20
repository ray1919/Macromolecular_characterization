library(shiny)
library(readr)
library(dplyr)
library(tagsinput)
library(DT)
library(openxlsx)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(ggplot2)

addResourcePath("www", "www")
ref <- list(
    HC = AAString("QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAIYPGNGDTSYNQKFKGKATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWGAGTTVTVSAASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKAEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"),
    LC = AAString("QIVLSQSPAILSASPGEKVTMTCRASSSVSYIHWFQQKPGSSPKPWIYATSNLASGVPVRFSGSGSGTSYSLTISRVEAEDAATYYCQQWTSNPPTFGGGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC")
)

ui <- fluidPage(
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'www/style.css'),
    navbarPage("大分子表征应用 v0.1",
               tabPanel(
                   "数据筛选",
                   sidebarLayout(
                       sidebarPanel(
                           sliderInput("confidence",
                                       "Confidence Score:",
                                       min = 1,
                                       max = 100,
                                       value = 80),
                           numericInput("delta",
                                        "Delta (ppm):",
                                        value = 5),
                           tagsTextInput("modification",
                                         "Modification:",
                                         value = "nonspecific, specific,None, Gln->Pyro-Glu, Lys_loss, G0F, G1F, G2F"),
                           fileInput("raw_csv",
                                     "CSV输入文件",
                                     accept = c(".csv"),
                                     buttonLabel = "选择"),
                           downloadButton("downloadData", "下载"),
                           tags$hr(),
                           tags$a(href="www/080_try.csv", target="_blank", "测试数据"),
                           tags$br(),
                           tags$img(src="BPSP_200x120.png", alt="生物医药公共服务平台", width="100%"),
                           width = 3
                       ),
                       mainPanel(
                           DT::dataTableOutput("flt_tbl"),
                           tags$h2("默认筛选逻辑："),
                           tags$ol(
                               tags$li("可信度大于等于80%；"),
                               tags$li("偏差小于等于5ppm；"),
                               tags$li("修饰字段包含nonspecific、specific、None、Gln->Pyro-Glu、Lys_loss、GOF、G1F、G2F；"),
                               tags$li("肽段字段去除重复。")
                           ),
                           width = 9
                           )
                   )
               ),
               tabPanel(
                   "覆盖率信息",
                   fluidRow(
                       
                       column(3,
                              tableOutput("aln"),
                              tags$hr(),
                              tableOutput("cov")
                       ),
                       
                       column(9,
                              plotOutput("covPlot")
                       )
                   )
               )
    )

)

server <- function(input, output) {

    raw_tbl <- reactive({
        rt <- tibble()
        if (!is.null(input$raw_csv)) {
            raw <- read_csv(input$raw_csv$datapath, show_col_types = F)
            all_mods <- unique(raw$Modification)
            mods <- str_split(input$modification, ", *") %>% unlist()
            pat_mods <- all_mods[grepl(paste(mods, collapse="|"), all_mods)]
            rt <- raw %>% 
                filter(`Confidence Score` >= input$confidence,
                       `Delta (ppm)` <= input$delta,
                       Modification %in% pat_mods) %>%
                distinct(`Peptide Sequence`, Protein, .keep_all = TRUE) %>%
                dplyr::select(`Peptide Sequence`,
                              `Modification`,
                              `M/Z`,
                              `Mono Mass Exp.`,
                              `Avg Mass Exp.`,
                              `Mono Mass Theo.`,
                              `RT`,
                              `Delta (ppm)`,
                              Protein)
        }
        rt
    })
    
    aln_tbl <- reactive({
        if (nrow(raw_tbl()) == 0)
            return(tibble(seqnames = character(),
                          start = character(),
                          end = character(),
                          width = character()))
        range_lst <- list()
        withProgress(message = '开始比对数据...',
                     detail = '比对肽段序列：',
                     value = 0,
                     max = 1, {
                         for ( i in 1:nrow(raw_tbl())) {
                            mp <- matchPattern(raw_tbl()$`Peptide Sequence`[i], ref[[raw_tbl()$Protein[i]]])
                            range_lst[[raw_tbl()$Protein[i]]][[i]] <- GRanges(seqnames = raw_tbl()$Protein[i], ranges = mp@ranges)
                        
                            incProgress(i/nrow(raw_tbl()), i)
                         }
        })
            
        ref_rng <- GRanges(seqnames = names(ref),
                           ranges = IRanges(1, width = sapply(ref, nchar)))
        
        b <- GRangesList(unlist(range_lst)) %>% unlist() %>% reduce()
        ss <- subsetByOverlaps(b, ref_rng)
        width(ss)
        
        as_tibble(ss) %>%
            mutate(ymin = if_else(seqnames == "HC", 1.5, 0.5),
                   ymax = if_else(seqnames == "HC", 2, 1))
    })
    
    cov_tbl <- reactive({
        if (nrow(aln_tbl()) == 0)
            return(tibble())
        aln_tbl() %>%
            group_by(seqnames) %>%
            summarise(length = sum(width)) %>%
            mutate(total = sapply(ref, nchar)[`seqnames`],
                   coverage = length / sapply(ref, nchar)[seqnames] * 100)
    })
    
    cov_plot <- reactive({
        if (nrow(aln_tbl()) == 0)
            return(ggplot())
        ggplot(aln_tbl()) +
            xlim(c(0,nchar(ref$HC))) +
            ylim(c(0,2)) + xlab("") + ylab("") +
            geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = seqnames), alpha = 0.7) +
            geom_rect(aes(xmin = 1, ymin = 1, xmax = nchar(ref$HC), ymax = 1.5), fill = "red3", alpha = 0.07) +
            geom_rect(aes(xmin = 1, ymin = 0, xmax = nchar(ref$LC), ymax = 0.5), fill = "deepskyblue", alpha = 0.07) +
            geom_text(aes(x = 225, y = 1.5, label = round(cov_tbl()$coverage[1],1))) +
            geom_text(aes(x = 106, y = 0.5, label = round(cov_tbl()$coverage[2],1))) +
            theme_bw() +
            theme(panel.border = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y = element_blank(),
                  legend.position = c(0.85, 0.3))
    })
    
    output$aln <- renderTable(aln_tbl() %>% dplyr::select(1:4))
    output$cov <- renderTable(cov_tbl())
    output$covPlot <- renderPlot(cov_plot())
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".xlsx", sep="")
        },
        content = function(file) {
            wb <- createWorkbook(creator = "Ryan Zhao")
            addWorksheet(wb, "Table1")
            writeDataTable(wb, 1, raw_tbl(), withFilter = F)
            saveWorkbook(wb, file, overwrite = T)
        }
    )
    
    output$flt_tbl <- renderDT(raw_tbl(), filter = 'top', options =
                                   list(pageLength = 15, scrollX = TRUE), 
                               rownames= FALSE)
    
}

# Run the application 
shinyApp(ui = ui, server = server)
