
# data folder
dir_home<- fs::path(here::here() %>% 
                      dirname())
dir_analysis <- dir_home |> dirname() |> dirname()

dir_data <- fs::path(dir_analysis,
                     "1_analysis_ready_data")

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
