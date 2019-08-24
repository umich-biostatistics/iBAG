
# Introduction page to be launched on startup of iBAG to guide the user
# to the various iBAG functions

introduction_page <- tabItem(
  
  # Name of this tab
  tabName = "introduction",
  
  # #HTML formatting adds gap between top of page
  # br(),
  # br(),
  
  # Markdown document of Introductory contents
  #includeMarkdown("iBAG/markdown/introduction_markdown.Rmd")
  includeMarkdown("markdown/introduction_markdown.Rmd"),
  
  # Image below the Introduction Markdown content
  column(
    width = 12, 
    offset = 1,
    shiny::img(src = "iBAG_Modeling_Strategy.png",
               width = 850,
               height = 370)
  )
  
)

