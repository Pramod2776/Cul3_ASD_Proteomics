setwd("/Users/pramod/Desktop/Proteomics/Fwd__SL_combat_script_validate/GO_enrichment_up_down/")
library(tidyverse)
options(tibble.width = Inf)
terms <- lapply(
  setNames(
    nm = readxl::excel_sheets("P7_CX.xlsx")
  ),
  function(sheet) {
    readxl::read_xlsx(
      "P7_CX.xlsx",
      sheet = sheet
    ) %>%
      mutate(Sheet = sheet)
  }
) %>%
  bind_rows()

go_plots<-lapply(
  setNames(
    nm = sort(unique(terms$Sheet))
  ),
  function(sheet){
    ggplot(
      data=terms %>%
        filter(sheet == sheet),
      mapping = aes(
        y = -log10(p.value) * Direction,
        x = reorder(
          str_wrap(term.name, 40),
          -log10(p.value) * Direction
        ),
        fill = as.character(Direction)
      )
    )+
      coord_flip()+
      geom_bar(
        stat = "identity", position = "dodge", colour = "black",
        show.legend = FALSE
      )+
      geom_hline(
        yintercept = c(-log10(0.05), log10(0.05)),
        colour = "red", linetype = "dashed"
      ) +
      labs(
        title = sheet,
        y = expression("-log"["10"]*"FDR")
      ) +
      scale_fill_manual(
        values = c(
          "-1" = "blue",
          "1" = "gold"
        )
      )+
      theme_bw() +
      theme(
        text = element_text(size = 30),
        axis.title.y = element_blank()
      )
  }
  
  
  )
pdf("DEG_GOEnrichmment_Selected_P7_CX.pdf", width = 16, height = 12)
go_plots
dev.off()





















