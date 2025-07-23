library(dplyr)

cia_wf_data <- read.csv("countries.csv")

colnames_for_anlysis <- c("Country","People.and.Society..Life.expectancy.at.birth...male",
                          "People.and.Society..Life.expectancy.at.birth...female",
                          "Economy..Real.GDP..purchasing.power.parity.")

cia_wf_data_le_vs_gdp <- cia_wf_data %>% select(all_of(colnames_for_anlysis))

index_2023 <- (grepl("(2023 est.)", cia_wf_data_le_vs_gdp$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))


cia_wf_data_2023 <- cia_wf_data_le_vs_gdp[index_2023,]

index_bil <- (grepl("billion (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))
index_tril <- (grepl("trillion (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))
index_mil <- (grepl("million (2023 est.)", cia_wf_data_2023$Economy..Real.GDP..purchasing.power.parity., fixed = TRUE))


colnames(cia_wf_data_2023) <- c("Country",
                              "Life_expectancy_M",
                              "Life_expectancy_F",
                              "GDP_PPP")

library(readr)

cia_wf_data_2023 = cia_wf_data_2023 %>%
  mutate(across(-c(Country),.fns = parse_number))

# commented as repating this will be problematic

# cia_wf_data_2023$GDP_PPP[index_bil] = 1000 * cia_wf_data_2023$GDP_PPP[index_bil]
# cia_wf_data_2023$GDP_PPP[index_tril] = 1000000 * cia_wf_data_2023$GDP_PPP[index_tril]

# cia_wf_data_2023 <- na.omit(cia_wf_data_2023)

plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Life_expectancy_M)
plot(log(cia_wf_data_2023$GDP_PPP),cia_wf_data_2023$Life_expectancy_F)

U1 = ecdf(cia_wf_data_2023$Life_expectancy_F)(cia_wf_data_2023$Life_expectancy_F)
U2 = ecdf(cia_wf_data_2023$Life_expectancy_M)(cia_wf_data_2023$Life_expectancy_M)

plot(U1,U2)
