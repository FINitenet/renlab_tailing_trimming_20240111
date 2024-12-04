# for tailing and trimming bubble plot
options(warn=-1)
library(ggplot2)
library(ggforce)
library(tidyverse)
library(patchwork)

files <- list.files(path = "/bios-store1/chenyc/Project/YHR/Project_yanghuiru_240510N/4_163.results/1D2_865_240510N_R1", pattern = "*.txt")

dir.create("5_GMC_analysis/plot_bubble")

buble_plot <- function(matrix.path, title, yanse) {
  df <- read.table(matrix.path, skip = 1)
  #df <- read.table("4_163.results/67_OXAP_R1_11x11/ath-MIR156a-3p_star.txt", skip = 1)
  sum <- sum(df)
  df <- sqrt(df / sum(df))*0.75 
  names(df) <- c(10:0)
  y0 <- c(10:0)
  df <- cbind(y0, df)
  df_long <- gather(df, x0, R, "10":"0")
  df_long$x0 <- as.integer(df_long$x0)

  ggplot(df_long) +
    geom_circle(aes(x0 = x0, y0 = y0, r = R), show.legend = FALSE, colour = NA, fill = yanse, na.rm = TRUE) +
    scale_x_reverse(limits = c(11, -1), breaks = seq(0, 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1, 11), breaks = seq(0, 10), expand = c(0, 0)) +
    # scale_size_area()+
    coord_fixed() +
    theme_bw() +
    theme(
      panel.grid.minor = element_line(colour = NA),
      panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12, face = "bold", colour = "black"),
      plot.title = element_text(size = 12, hjust = 0.5, face = "italic")
    ) +
    xlab("Length of trimming") +
    ylab("Length of tailing") +
    ggtitle(paste0(title, "(", sum, ")"))
}

buble_plot1 <- function(matrix.path, title, yanse) {
  df <- read.table(matrix.path, skip = 1)
  #df <- read.table("4_163.results/67_OXAP_R1_11x11/ath-MIR156a-3p_star.txt", skip = 1)
  sum <- sum(df)
  df <- sqrt(df / sum(df))*0.75 
  names(df) <- c(10:0)
  y0 <- c(10:0)
  df <- cbind(y0, df)
  df_long <- gather(df, x0, R, "10":"0")
  df_long$x0 <- as.integer(df_long$x0)
  
  ggplot(df_long) +
    geom_circle(aes(x0 = x0, y0 = y0, r = R), show.legend = FALSE, colour = NA, fill = yanse, na.rm = TRUE) +
    scale_x_reverse(limits = c(11, -1), breaks = seq(0, 10), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1, 11), breaks = seq(0, 10), expand = c(0, 0)) +
    # scale_size_area()+
    coord_fixed() +
    theme_bw() +
    theme(
      panel.grid.minor = element_line(colour = NA),
      panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
      panel.background = element_rect(fill = NA),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 9, face = "bold"),
      axis.text = element_text(size = 9, face = "bold", colour = "black"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7, face = "bold", hjust = 0.5),
      plot.background = element_rect(colour = NA),
      legend.position = "none") +
    xlab("Length of trimming") +
    ylab("Length of tailing") +
    ggtitle(title, paste0("total reads: ",sum))
}


for (i in seq_along(files)) { 
  miRNA_title <- gsub("\\.txt", "", files[i])
  miRNA_id <- gsub("_star", "\\*", miRNA_title)
  print(paste0("Processing(", i, "/", length(files), "):", miRNA_id))
  P1 <- buble_plot1(matrix.path = paste0("4_163.results/hen1_8_240510N_R1/", files[i]), title = "hen1_8_240510N_R1", yanse = "#9b4b8f")
  P2 <- buble_plot1(matrix.path = paste0("4_163.results/1D2_865_240510N_R1/", files[i]), title = "1D2_865_240510N_R1", yanse = "#9b4b8f")
  P3 <- buble_plot1(matrix.path = paste0("4_163.results/1F2_576_240510N_R1/", files[i]), title = "1F2_576_240510N_R1", yanse = "#9b4b8f")
  P4 <- buble_plot1(matrix.path = paste0("4_163.results/1P1_867_240510N_R1/", files[i]), title = "1P1_867_240510N_R1", yanse = "#9b4b8f")
  P5 <- buble_plot1(matrix.path = paste0("4_163.results/2F_4_4138_240510N_R1/", files[i]), title = "2F_4_4138_240510N_R1", yanse = "#9b4b8f")
  P6 <- buble_plot1(matrix.path = paste0("4_163.results/2L_14_15_496_240510N_R1/", files[i]), title = "2L_14_15_496_240510N_R1", yanse = "#9b4b8f")
  P7 <- buble_plot1(matrix.path = paste0("4_163.results/2L_896_240510N_R1/", files[i]), title = "2L_896_240510N_R1", yanse = "#9b4b8f")
  merge <-(P1 | ( (P2 + P3 + P4) / (P5 + P6 + P7)))  + plot_annotation(title = miRNA_id, theme = theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold")))
  ggsave(paste0("5_GMC_analysis/plot_bubble/", miRNA_title, ".pdf"), merge, width = 12, height = 6, units = "in", dpi = 300)
}
