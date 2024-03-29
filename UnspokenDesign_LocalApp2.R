# libraries
library(shiny)
library(tidyverse)
library(dplyr)
library(reshape2)
library(shinycssloaders)
library(shinythemes)
library(DT)
library(tidyr)

source("https://raw.githubusercontent.com/klattery/Unspoken/master/Attraction_Conversion_Design_Generator_20SET20_kl3.R")
source("https://raw.githubusercontent.com/klattery/Unspoken/master/UnspokenDesign_RShiny2.R")
shinyApp(ui,server)
