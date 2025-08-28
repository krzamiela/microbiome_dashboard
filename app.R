library(tidyverse)
library(shiny)
library(bslib)
library(gt)
library(patchwork)
library(scales)

########################
## Analysis Functions ##
########################

# Group samples by binary cancer type and genus
# only cancer or no cancer (BC/MBC1 or BC/MBC2)
get_binary_groups <- function(region_data) {
  region_data %>%
    pivot_longer(cols = matches("^(BC|MBC)"),
                 names_to = "Sample",
                 values_to = "Count") %>%
    select(Sample, Taxa, Kingdom, Phylum, Class, Order, Family, Genus, Species, Count) %>%
    filter(!str_detect(Sample, "^MBC[3-5]")) %>%
    mutate(Taxon = paste(Genus, Species, sep = " ")) %>%
    mutate(cancer_status = case_when(
      str_detect(Sample, "^BC1")  ~ "breast cancer",
      str_detect(Sample, "^MBC1") ~ "breast cancer",
      str_detect(Sample, "^BC2")  ~ "no breast cancer",
      str_detect(Sample, "^MBC2") ~ "no breast cancer",
      TRUE ~ "unknown"))
}

# Group samples by all cancer types and genus
get_all_groups <- function(region_data) {
  region_data %>%
    pivot_longer(cols = matches("^(BC|MBC)"),
                 names_to = "Sample",
                 values_to = "Count") %>%
    select(Sample, Taxa, Kingdom, Phylum, Class, Order, Family, Genus, Species, Count) %>%
    mutate(Taxon = paste(Genus, Species, sep = " ")) %>%
    mutate(cancer_status = case_when(
      str_detect(Sample, "^BC1")  ~ "BC",
      str_detect(Sample, "^MBC1") ~ "BC",
      str_detect(Sample, "^BC2")  ~ "No BC",
      str_detect(Sample, "^MBC2") ~ "No BC",
      str_detect(Sample, "^MBC3") ~ "No BC, abnormal mammogram",
      str_detect(Sample, "^MBC4") ~ "No BC, 1st deg BC",
      str_detect(Sample, "^MBC5") ~ "No BC, 1st deg BC in study",
      TRUE ~ "unknown")) %>%
    filter(cancer_status != "No BC, 1st deg BC in study") # sample size too small, skews the data
}

# get unclassified organisms at all taxa levels
get_unclassified_all <- function(qiime_bc_taxa) {
  qiime_bc_taxa %>%
    filter(Kingdom == "Unclassified" | Phylum == "Unclassified" | Class == "Unclassified" | Order == "Unclassified" | Family == "Unclassified" | Genus == "Unclassified" | Species == "Unclassified") %>%
    mutate(unclassified_category = case_when(
      Kingdom == "Unclassified" ~ "Kingdom",
      Phylum == "Unclassified" ~ "Phylum",
      Class == "Unclassified" ~ "Class",
      Order == "Unclassified" ~ "Order",
      Family == "Unclassified" ~ "Family",
      Genus == "Unclassified" ~ "Genus",
      Species == "Unclassified" ~ "Species",
      TRUE ~ "Unknown"
    )) %>%
    group_by(unclassified_category) %>%
    summarize(unclass_count = n(),
              unclass_perc = n() / nrow(qiime_bc_taxa),
              unclass_org_count = sum(`Total Count`)) %>%
    ungroup()
}

# calculate prevalence at the genus level
# expects taxon type to be capitalized (Species, Genus, etc)
prevalence <- function(taxon_name, taxon_type, sample_data) {
  if (taxon_type == "Species") {
    result <- sample_data %>%
      filter(Taxon == taxon_name) %>%
      group_by(Sample, cancer_status) %>%
      summarize(present = as.integer(any(Count > 0))) %>%
      ungroup() %>%
      group_by(cancer_status) %>%
      summarize(
        prevalence = sum(present) / n(),
        total_samples = n()) %>%
      ungroup() %>%
      mutate(y_label = str_c(round(prevalence * 100, 2), "%"))
  }  else {
    result <- sample_data %>%
      filter(.data[[taxon_type]] == taxon_name) %>%
      group_by(Sample, cancer_status) %>%
      summarize(present = as.integer(any(Count > 0))) %>%
      ungroup() %>%
      group_by(cancer_status) %>%
      summarize(
        prevalence = sum(present) / n(),
        total_samples = n()) %>%
      ungroup() %>%
      mutate(y_label = str_c(round(prevalence * 100, 2), "%"))
  }
  
  return(result)
}

# calculate relative abundance
rel_abundance <- function(taxon_name, taxon_type, sample_data) {
  if (taxon_type == "Species") {
    result <- sample_data %>%
      group_by(Sample, cancer_status) %>%
      summarize(abundance_tax = sum(Count[Taxon == taxon_name]),
                total_abundance_sample = sum(Count),
                rel_abundance_tax = sum(Count[Taxon == taxon_name]) / sum(Count)) %>%
      ungroup()
  } else {
    result <- sample_data %>%
      group_by(Sample, cancer_status) %>%
      summarize(abundance_tax = sum(Count[.data[[taxon_type]] == taxon_name]),
                total_abundance_sample = sum(Count),
                rel_abundance_tax = sum(Count[.data[[taxon_type]] == taxon_name]) / sum(Count)) %>%
      ungroup()
  }
  
  return(result)
}

###########################
## Define plot functions ##
###########################

# prevalence plot
prev_plot <- function(plot_data, taxon, flipxy, region) {
  
  # create plot
  if (flipxy == "no") {
    plot <- ggplot(data = plot_data, aes(x = cancer_status, y = prevalence, fill = cancer_status)) + 
      geom_col(alpha = 0.8, 
               show.legend = FALSE) +
      scale_y_continuous(name = "Prevalence", 
                         limits = c(0, 1), 
                         labels = label_percent(scale = 100)) +
      scale_fill_brewer(palette = "Dark2") +
      geom_text(aes(label = y_label), 
                vjust = -1, 
                size = 3.5, 
                color = "black") +
      labs(x = NULL) +
      theme_minimal() +
      ggtitle(str_c(region, ": Prevalence of ", taxon)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(plot)
  }
  else {
    plot <- ggplot(data = plot_data, aes(x = cancer_status, y = prevalence, fill = cancer_status)) + 
      geom_col(alpha = 0.8, show.legend = FALSE) +
      scale_y_continuous(name = "Prevalence", 
                         limits = c(0, 1), 
                         labels = label_percent(scale = 100)) +
      scale_fill_brewer(palette = "Dark2") +
      geom_text(aes(label = y_label), 
                hjust = 1.25, 
                size = 3.5, 
                color = "black") +
      coord_flip() +
      labs(x = NULL) +
      theme_minimal() +
      ggtitle(str_c(region, ": Prevalence of ", taxon)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(plot)
  }
}

# abundance plot
ab_plot <- function(plot_data, taxon, flipxy, region) {
  
  # create plot
  if (flipxy == "no") {
    plot <- ggplot(data = plot_data, mapping = aes(x = cancer_status, y = rel_abundance_tax, fill = cancer_status)) + 
      geom_boxplot(alpha = 0.7,
                   outlier.shape = NA,
                   show.legend = FALSE) +
      geom_jitter(alpha = 0.6,
                  color = "black",
                  width = 0.2,
                  show.legend = FALSE) +
      scale_y_continuous(name = "Relative Abundance",
                         labels = label_percent()) +
      scale_fill_brewer(palette = "Dark2") +
      labs(x = NULL) +
      theme_minimal() +
      ggtitle(str_c(region, ": Relative ", taxon, " Abundance")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(plot)
  }
  else {
    plot <- ggplot(data = plot_data, mapping = aes(x = cancer_status, y = rel_abundance_tax, fill = cancer_status)) + 
      geom_boxplot(alpha = 0.7,
                   outlier.shape = NA,
                   show.legend = FALSE) +
      geom_jitter(alpha = 0.6,
                  color = "black",
                  width = 0.2,
                  show.legend = FALSE) +
      scale_y_continuous(name = "Relative Abundance",
                         labels = label_percent()) +
      scale_fill_brewer(palette = "Dark2") +
      coord_flip() +
      labs(x = NULL) +
      theme_minimal() +
      ggtitle(str_c(region, ": Relative ", taxon, " Abundance")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    return(plot)
  }
}


# Define UI for application that draws plots
ui <- tagList(
  accordion(
    id = "instructions_accordion",
    accordion_panel(
      title = "How to use this dashboard",
      HTML(
        "<ul>
              <li>First, upload your V1-V3 and V3-V4 CSV files using the panels on the left.</li>
              <li>Next, select a taxonomic rank and specific taxon to analyze.</li>
              <li>Choose whether to group samples into <em>binary groups</em> (cancer vs. no cancer) or view <em>all groups</em>.</li>
              <li>Note: Some selections may take time to render due to large dataset sizes.</li>
            </ul>
            <p style='margin-top: 1rem; font-size: 0.9em; color: #555;'>
              <em>Data source: Dr. Ece Mutlu’s Lab, University of Illinois Chicago</em>
            </p>"
      )
    )
  ),
  page_sidebar(
    title = "Breast Cancer 16S Microbiome Data",
    sidebar = sidebar(
      card(
        card_header("Upload V1-V3 file"),
        fileInput("v1_v3_file",
                  label = NULL)
      ),
      card(
        card_header("Upload V3-V4 file"),
        fileInput("v3_v4_file",
                  label = NULL)
      ),
      card(
        card_header("Taxonomic Rank"),
        selectInput(
          inputId = "tax_select",
          label = NULL,
          choices = c("Phylum" = "Phylum",
                      "Class" = "Class",
                      "Order" = "Order",
                      "Family" = "Family",
                      "Genus" = "Genus"),
          selectize = FALSE
        )
      ),
      card(
        card_header("Entity"),
        uiOutput("entity_dropdown")
      ),
      card(
        card_header("Analysis Type"),
        radioButtons(
          inputId = "radio_analysis", 
          label = NULL, 
          choices = c( "Binary Groups" = 1, 
                      "All Groups" = 2)
        )
      )
    ),
    card(
      card_header("Global Metrics"),
      layout_column_wrap(
        width = 2,
        style = "gap: 1rem;",
        card(
          plotOutput("globePlot", height = "900px", width = "100%")
        ),
        card(
          gt_output("uc_tbl")
        )
      )
    ),
    navset_card_underline(
      title = "Analysis",
      nav_panel("Prevalence", plotOutput("prevPlot", height = "900px")),
      nav_panel("Abundance", plotOutput("abPlot", height = "900px"))
    )
  )
)

# Define server logic for plotting functions
server <- function(input, output) {
  options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB limit, adjust as needed

  # Read files from uploaded paths
  v1_v3_data <- reactive({
    req(input$v1_v3_file)
    read_csv(input$v1_v3_file$datapath)
  })
  
  v3_v4_data <- reactive({
    req(input$v3_v4_file)
    read_csv(input$v3_v4_file$datapath)
  })
  
  # get unclassified organisms
  unclass_all_data <- reactive({
    req(v1_v3_data(), v3_v4_data())
    
    unclass_v1_v3 <- get_unclassified_all(v1_v3_data())
    unclass_v3_v4 <- get_unclassified_all(v3_v4_data())
    
    # merge into one df
    unclass_all <- merge(unclass_v1_v3, unclass_v3_v4, by = "unclassified_category", suffixes = c("_v1v3", "_v3v4")) %>%
      select(unclassified_category, unclass_perc_v1v3, unclass_perc_v3v4) %>%
      pivot_longer(cols = c(unclass_perc_v1v3, unclass_perc_v3v4), names_to = "region", values_to = "unclass_perc") %>%
      mutate(region = recode(region,
                             unclass_perc_v1v3 = "V1-V3",
                             unclass_perc_v3v4 = "V3-V4"),
             unclassified_category = factor(unclassified_category, levels = rev(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))))
    
    unclass_all
  })
  
  # get list of taxa
  filtered_taxa <- reactive({
    req(v1_v3_data(), v3_v4_data())

    v1_v3_taxa <- v1_v3_data() %>%
      select(all_of(input$tax_select)) %>%
      unique()
    
    v3_v4_taxa <- v3_v4_data() %>%
      select(all_of(input$tax_select)) %>%
      unique()
    
    # merge and only keep unique
    unique_taxa <- bind_rows(v1_v3_taxa, v3_v4_taxa) %>%
      distinct() %>%
      arrange(across(everything())) %>%
      pull() %>%
      str_trim() %>%
      str_replace_all('^"|"$', '')
    
    # create list
    setNames(unique_taxa, unique_taxa)
  })
  
  # render taxa dropdown
  output$entity_dropdown <- renderUI({
    req(filtered_taxa())
    
    selectInput(
      inputId = "entity_select",
      label = NULL,
      choices = filtered_taxa(),
      selectize = FALSE
    )
  })
  
  # binary groups for V1-V3 and V3-V4
  reg_1_sample_taxa_bin <- reactive({
    req(v1_v3_data())
    get_binary_groups(v1_v3_data())
  })
  
  reg_2_sample_taxa_bin <- reactive({
    req(v3_v4_data())
    get_binary_groups(v3_v4_data())
  })
  
  # all groups for V1-V3 and V3-V4
  reg_1_sample_taxa_group <- reactive({
    req(v1_v3_data())
    get_all_groups(v1_v3_data())
  })
  
  reg_2_sample_taxa_group <- reactive({
    req(v3_v4_data())
    get_all_groups(v3_v4_data())
  })
  
  # global unclassified
  output$globePlot <- renderPlot({
    req(unclass_all_data())
    
    # plot
    uc_plot <- ggplot(data = unclass_all_data(), mapping = aes(x = unclassified_category, y = unclass_perc, fill = region)) +
      geom_col(position = "dodge",
               alpha = 0.8) +
      geom_text(mapping = aes(label = percent(unclass_perc, accuracy = 0.1)),
                position = position_dodge(width = 0.9),
                vjust = 0.5,
                hjust = -0.1,
                size = 3) +
      scale_y_continuous(limits = c(0, 1),
                         labels = scales::percent_format(accuracy = 1)) +
      scale_fill_brewer(palette = "Dark2") +
      coord_flip() +
      labs(title = "Percentage Unclassified",
           x = NULL,
           y = NULL,
           fill = "Sequencing Region") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "top")
    
    uc_plot
  })
  
  # global table
  output$uc_tbl <- render_gt({
    req(v1_v3_data(), v3_v4_data())
    
    unclass_v1_v3 <- get_unclassified_all(v1_v3_data())
    unclass_v3_v4 <- get_unclassified_all(v3_v4_data())

    # keep table of counts
    uc_tbl <- merge(unclass_v1_v3, unclass_v3_v4, by = "unclassified_category", suffixes = c("_v1v3", "_v3v4")) %>%
      select(-matches("_perc")) %>%
      mutate(unclassified_category = factor(
        unclassified_category,
        levels = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      )) %>%
      arrange(unclassified_category) %>%
      gt() %>%
      tab_header(title = "Organism Counts") %>%
      cols_label(
        unclassified_category = "Category",
        unclass_count_v1v3 = "Unclassified",
        unclass_org_count_v1v3 = "Total",
        unclass_count_v3v4 = "Unclassified",
        unclass_org_count_v3v4 = "Total"
      ) %>%
      tab_spanner(
        label = "V1-V3",
        columns = c(unclass_count_v1v3, unclass_org_count_v1v3)
      ) %>%
      tab_spanner(
        label = "V3-V4",
        columns = c(unclass_count_v3v4, unclass_org_count_v3v4)
      ) %>%
      tab_style(
        style = cell_fill(color = "#ccece6"),
        locations = cells_body(columns = c(
          unclass_count_v1v3,
          unclass_org_count_v1v3
        )) # v1-v3 style
      ) %>%
      tab_style(
        style = cell_fill(color = "#fee6ce"),
        locations = cells_body(columns = c(
          unclass_count_v3v4,
          unclass_org_count_v3v4
        )) # v3-v4 style
      )
    
    uc_tbl
  })

  output$prevPlot <- renderPlot({
    req(v1_v3_data(), v3_v4_data(), input$entity_select, input$tax_select, input$radio_analysis)

    if (input$radio_analysis == 1) {
      # binary group analysis
      # calculate prevalence
      reg_1_prev_bin <- prevalence(input$entity_select, input$tax_select, reg_1_sample_taxa_bin())
      
      # generate plot with significance
      reg_1_plot <- prev_plot(reg_1_prev_bin, input$entity_select, "no", "V1-V3")
      
      # calculate prevalence
      reg_2_prev_bin <- prevalence(input$entity_select, input$tax_select, reg_2_sample_taxa_bin())
      
      # generate plot with significance
      reg_2_plot <- prev_plot(reg_2_prev_bin, input$entity_select, "no", "V3-V4")

      reg_1_plot + reg_2_plot
      
    } else if (input$radio_analysis == 2) {
      # all groups
      # prevalence
      reg_1_prev_group <- prevalence(input$entity_select, input$tax_select, reg_1_sample_taxa_group())
      
      # generate plot
      reg_1_plot <- prev_plot(reg_1_prev_group, input$entity_select, "yes", "V1-V3")
      
      # prevalence
      reg_2_prev_group <- prevalence(input$entity_select, input$tax_select, reg_2_sample_taxa_group())
      
      # generate plot
      reg_2_plot <- prev_plot(reg_2_prev_group, input$entity_select, "yes", "V3-V4")
      
      # display plot
      reg_1_plot + reg_2_plot
      
    } else {
      # binary group analysis
      # calculate prevalence
      reg_1_prev_bin <- prevalence(input$entity_select, input$tax_select, reg_1_sample_taxa_bin())
      
      # generate plot with significance
      reg_1_plot <- prev_plot(reg_1_prev_bin, input$entity_select, "no", "V1-V3")
      
      # calculate prevalence
      reg_2_prev_bin <- prevalence(input$entity_select, input$tax_select, reg_2_sample_taxa_bin())
      
      # generate plot with significance
      reg_2_plot <- prev_plot(reg_2_prev_bin, input$entity_select, "no", "V3-V4")
      
      reg_1_plot + reg_2_plot
    }
  })
  
  output$abPlot <- renderPlot({
    req(v1_v3_data(), v3_v4_data(), input$entity_select, input$tax_select, input$radio_analysis)
    
    if (input$radio_analysis == 1) {
      
      # binary groups
      # V1-V3
      reg_1_abundance_bin <- rel_abundance(input$entity_select, input$tax_select, reg_1_sample_taxa_bin())
      
      # Remove the 0s
      reg_1_ab_bin_filtered <- reg_1_abundance_bin %>%
        filter(rel_abundance_tax > 0)

      # plot
      reg_1_plot <- ab_plot(reg_1_ab_bin_filtered, input$entity_select, "no", "V1-V3")
      
      # V3-V4
      reg_2_abundance_bin <- rel_abundance(input$entity_select, input$tax_select, reg_2_sample_taxa_bin())
      
      # Remove the 0s
      reg_2_ab_bin_filtered <- reg_2_abundance_bin %>%
        filter(rel_abundance_tax > 0)
      
      # plot
      reg_2_plot <- ab_plot(reg_2_ab_bin_filtered, input$entity_select, "no", "V3-V4")
      
      # display
      reg_1_plot + reg_2_plot
      
    } else if (input$radio_analysis == 2) {
      
      # all groups
      # V1-V3
      reg_1_abundance_group <- rel_abundance(input$entity_select, input$tax_select, reg_1_sample_taxa_group())
      
      # Remove the 0s
      reg_1_ab_grp_filtered <- reg_1_abundance_group %>%
        filter(rel_abundance_tax > 0)
      
      # plot
      reg_1_plot <- ab_plot(reg_1_ab_grp_filtered, input$entity_select, "yes", "V1-V3")
      
      # V3-V4
      reg_2_abundance_group <- rel_abundance(input$entity_select, input$tax_select, reg_2_sample_taxa_group())
      
      # Remove the 0s
      reg_2_ab_grp_filtered <- reg_2_abundance_group %>%
        filter(rel_abundance_tax > 0)
      
      # plot
      reg_2_plot <- ab_plot(reg_2_ab_grp_filtered, input$entity_select, "yes", "V3-V4")
      
      # display
      reg_1_plot + reg_2_plot

    } else {
      
      # binary groups
      # V1-V3
      reg_1_abundance_bin <- rel_abundance(input$entity_select, input$tax_select, reg_1_sample_taxa_bin())
      
      # Remove the 0s
      reg_1_ab_bin_filtered <- reg_1_abundance_bin %>%
        filter(rel_abundance_tax > 0)
      
      # plot
      reg_1_plot <- ab_plot(reg_1_ab_bin_filtered, input$entity_select, "no", "V1-V3")
      
      # V3-V4
      reg_2_abundance_bin <- rel_abundance(input$entity_select, input$tax_select, reg_2_sample_taxa_bin())
      
      # Remove the 0s
      reg_2_ab_bin_filtered <- reg_2_abundance_bin %>%
        filter(rel_abundance_tax > 0)
      
      # plot
      reg_2_plot <- ab_plot(reg_2_ab_bin_filtered, input$entity_select, "no", "V3-V4")
      
      # display
      reg_1_plot + reg_2_plot
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
