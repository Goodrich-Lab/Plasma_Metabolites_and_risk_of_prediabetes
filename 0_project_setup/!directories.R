

# data folder
dir_data <- fs::path("~",
                     "Desktop",
                     "progress",
                     "Data",
                     "0_Data_mirror_do_not_edit")

dir_home<- fs::path(here::here() %>% 
                          dirname())

# analysis ready data folder
dir_data_analysis <- fs::path(
  here::here() %>% 
  dirname(), "0_Data", "analysis ready data")

# Result folder
dir_results <- fs::path(
  here::here() %>%
    dirname(), "2_Results"
)

# Report folder
dir_reports <- fs::path(
  here::here() %>%
    dirname(), "3_Reports"
)
  
# CHS analysis ready data folder
dir_chs <- fs::path(
  here::here() %>%
    dirname() %>% 
    dirname() %>%
    dirname(), "1_analysis_ready_data"
    
)
