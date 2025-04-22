setwd("~/Dropbox/projects_WORKING/temporal_signal_comment/measurable_evolution/")
library(NELSI)

plot.results.cis <- function(trees, posterior_files, prior_files, labels = T, 
                             phylo_threshold=20, folder_dir = "log_files/", 
                             main_text = "", ...){
  if(labels){
    labels_text <- c("Time (years BP)", expression("Evol. rate"~10^{-5}~"subs/site/year"), 
                     "Error in tree height")
  }else{
    labels_text <- c("", "", "")
  }
  clock_CI_table <- matrix(NA, nrow = length(trees), ncol = 6)
  colnames(clock_CI_table) <- c("mean_posterior", "lower_posterior", "upper_posterior", 
                              "mean_prior", "lower_prior", "upper_prior")
  th_CI_table <- matrix(NA, nrow = length(trees), ncol = 6)
  colnames(th_CI_table) <- c("mean_posterior", "lower_posterior", "upper_posterior", 
                           "mean_prior", "lower_prior", "upper_prior")
  for(i in 1:length(trees)){
    temp_tree <- trees[[i]]
    true_tree_height <- max(allnode.times(temp_tree))
    simulation_replicate <- paste0("_sim", i, "_")
    print(simulation_replicate)
    temp_posterior <- tail(read.table(paste0(folder_dir, 
                                       grep(simulation_replicate, posterior_files, value = T)), 
                                 head = T), 1000)
    temp_prior <- tail(read.table(paste0(folder_dir, 
                                    grep(simulation_replicate, prior_files, value = T)), 
                              head = T), 1000)
    th_CI_table[i, 1:3] <- c(mean(true_tree_height - temp_posterior$tree.height), 
                            quantile(true_tree_height - temp_posterior$tree.height, c(0.025, 0.975))) / true_tree_height
    th_CI_table[i, 4:6] <- c(mean(true_tree_height - temp_prior$tree.height),
                            quantile(true_tree_height - temp_prior$tree.height, c(0.025, 0.975))) / true_tree_height
    clock_CI_table[i, 1:3] <- c(mean(temp_posterior$rate.mean), 
                                quantile(temp_posterior$rate.mean, c(0.025, 0.975)))
    clock_CI_table[i, 4:6] <- c(mean(temp_prior$rate.mean),
                                quantile(temp_prior$rate.mean, c(0.025, 0.975)))
  }

  clock_CI_table <- clock_CI_table[order(clock_CI_table[, 1]), ]
  par(mar = c(0.1, 5, 5, 2))
  plot(0, 0, ylim = c(-1, 4.3), xlim = c(0, 100), bty = "n", xaxt = "n", yaxt = "n", 
       xlab = "", type = "n", ylab = labels_text[1], main = main_text, ...)
  sampling_height <- max(allnode.times(trees[[1]], tipsonly = T, reverse = F))
  axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
       labels = c("0", "20", "1K", "10K", "20K"), las = 2)
  plot.tree.lines(trees[[1]], line.type = "l", rotation.angle = 3*pi/2, log.scale = T)
  lines(c(-10, 100), rep(log10(sampling_height), 2), lty = 1, lwd = 1, col = "red")
  lines(c(-10, 100), rep(log10(20), 2), lty = 2, lwd = 1, col = "purple")
  
  par(mar = c(1, 5, 0.1, 3))
  plot(c(-1, 101), log10(c(5e-7, 1e-4)), type = "n", bty = "n", yaxt = "n", 
       xaxt = "n", xlab = "", ylab  = labels_text[2])
  axis(2, at = log10(c(5e-7, 2e-6, 5e-6, 1.5e-5, 5e-5, 1e-4)), 
       labels = c("0.05", "0.2", "0.5", "1.5", "5.0",  "100"), las = 2)
  axis(1, at = c(0, 101), labels = c("", ""))
  #
  for(i in 1:nrow(clock_CI_table)){
    lines(c(i, i)-0.1, log10(clock_CI_table[i, 5:6]), col = "orange", lwd = 0.7)
    points(i-0.1, log10(clock_CI_table[i, 4]), pch = 20, col = "orange", cex = 0.5)
    lines(c(i, i)+0.1, log10(clock_CI_table[i, 2:3]), col = "blue", lwd = 1.2)
    points(i+0.1, log10(clock_CI_table[i, 1]), pch = 20, col = "blue")
  }
  lines(c(-5, 105), rep(log10(1.5e-5), 2), lty = 2, lwd = 2)

  th_CI_table <- th_CI_table[order(th_CI_table[, 1]), ]
  plot(c(-1, 101), c(-10, 2), type = "n", bty = "n", xaxt = "n", xlab = "", 
     yaxt = "n", ylab = labels_text[3])
  axis(1, at = c(0, 101), labels = c("", ""))
  axis(2, at = c(-10, -5, -2, 0, 2), labels = c("-10", "-5", "-2", "0", "2"), las = 2)
  for(i in 1:nrow(th_CI_table)){
    lines(c(i, i)-0.1, th_CI_table[i, 5:6], col = "orange", lwd = 0.7)
    points(i-0.1, th_CI_table[i, 4], pch = 20, col = "orange", cex = 0.5)
    lines(c(i, i)+0.1, th_CI_table[i, 2:3], col = "blue", lwd = 1.2)
    points(i+0.1, th_CI_table[i, 1], pch = 20, col = "blue")
  }
  lines(c(-5, 105), rep(0, 2), lty = 2, lwd = 2)
}

###########################################
########################################### 
#Plot with reasonable prior
###########################################
########################################### 

pdf("summary_all_estimates_correct_prior.pdf", width = 11, height = 7)
par(mfcol = c(3, 5))
#
trees <- read.tree("coalescent_late_sampling_dated_ultrametric.trees")
posterior_files <- dir(path = "log_files/", pattern = ".+ultra.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "log_files/", pattern = ".+ultra.+sample_prior.log")
par(mar = c(0.1, 5, 5, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = T, 
  main_text = expression("(a) Sampling window=0\n(ultrametric)"), cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_dated_0.5X_pd_threshold.trees")
posterior_files <- dir(path = "log_files/", pattern = ".+0.5X.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "log_files/", pattern = ".+0.5X.+sample_prior.log")
par(mar = c(0.1, 5, 5, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(b) Sampling window=\n0.5*phylo. threshold"), cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_dated.trees")
posterior_files <- dir(path = "log_files/", pattern = ".+under.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "log_files/", pattern = ".+under.+sample_prior.log")
par(mar = c(0.1, 5, 5, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(c) Sampling window=\nphylo. threshold"), cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_2X_pd_threshold_dated.trees")
posterior_files <- dir(path = "log_files/", pattern = ".+2X.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "log_files/", pattern = ".+2X.+sample_prior.log")
par(mar = c(0.1, 5, 5, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F,
  main_text = expression("(d) Sampling window=\n10*phylo. threshold"), cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_10X_pd_threshold_dated.trees")
posterior_files <- dir(path = "log_files/", pattern = ".+10X.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "log_files/", pattern = ".+10X.+sample_prior.log")
par(mar = c(0.1, 5, 5, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(e) Sampling window=\n100 * phylo. threshold"), cex.main = 1.2)
dev.off()

## Need to add main labels below!
###########################################
########################################### 
#Plot with misleading prior
###########################################
########################################### 

pdf("summary_all_estimates_misleading_prior.pdf", width = 11, height = 7)
par(mfcol = c(3, 5))
#
trees <- read.tree("coalescent_late_sampling_dated_ultrametric.trees")
posterior_files <- dir(path = "misleading_prior/", pattern = ".+ultra.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "misleading_prior/", pattern = ".+ultra.+sample_prior.+log")
par(mar = c(0.1, 5, 0.1, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = T, , 
  main_text = expression("(a) Sampling window=0\n(ultrametric)"),
  folder_dir = "misleading_prior/", cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_dated_0.5X_pd_threshold.trees")
posterior_files <- dir(path = "misleading_prior/", pattern = ".+0.5X.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "misleading_prior/", pattern = ".+0.5X.+sample_prior.+log")
par(mar = c(0.1, 5, 0.1, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(b) Sampling window=\n0.5*phylo. threshold"),
  folder_dir = "misleading_prior/", cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_dated.trees")
posterior_files <- dir(path = "misleading_prior/", pattern = ".+under.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "misleading_prior/", pattern = ".+under.+sample_prior.+log")
par(mar = c(0.1, 5, 0.1, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(c) Sampling window=\nphylo. threshold"),
  folder_dir = "misleading_prior/", cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_2X_pd_threshold_dated.trees")
posterior_files <- dir(path = "misleading_prior/", pattern = ".+2X.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "misleading_prior/", pattern = ".+2X.+sample_prior.+log")
par(mar = c(0.1, 5, 0.1, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(d) Sampling window=\n10*phylo. threshold"),
  folder_dir = "misleading_prior/", cex.main = 1.2)
#
trees <- read.tree("coalescent_late_sampling_10X_pd_threshold_dated.trees")
posterior_files <- dir(path = "misleading_prior/", pattern = ".+10X.+log")
posterior_files <- posterior_files[!grepl("sample_prior", posterior_files)]
prior_files <- dir(path = "misleading_prior/", pattern = ".+10X.+sample_prior.+log")
par(mar = c(0.1, 5, 0.1, 1))
plot.results.cis(trees, posterior_files, prior_files, labels = F, 
  main_text = expression("(e) Sampling window=\n100 * phylo. threshold"),
  folder_dir = "misleading_prior/", cex.main = 1.2)
dev.off()

