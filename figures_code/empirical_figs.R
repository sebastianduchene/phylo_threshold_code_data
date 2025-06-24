library(NELSI)
setwd("/home/sebastiand/Dropbox/projects_WORKING/temporal_signal_comment/HBV/analyses_julia")


##############################################
###### Phylodyn threshold

phylo_threshold_hcc_trees <- c("T0.trees.hcc.tre", "T500.trees.hcc.tre", 
                        "T1000.trees.hcc.tre", "T_3020.trees.hcc.tre")

T0_hbv_prior <- tail(read.table("T0.log", head = T), 1000)
T0_hbv_sample_prior <- tail(read.table("T0_sample_prior.log", head = T), 1000)

T500_hbv_prior <- tail(read.table("T500.log", head = T), 1000)
T500_hbv_sample_prior <- tail(read.table("T500_sample_prior.log", head = T), 1000)

T1000_hbv_prior <- tail(read.table("T1000.log", head = T), 1000)
T1000_hbv_sample_prior <- tail(read.table("T1000_sample_prior.log", head = T), 1000)

T_3020_hbv_prior <- tail(read.table("T_3020.log", head = T), 1000)
T_3020_hbv_sample_prior <- tail(read.table("T_3020_sample_prior.log", head = T), 1000)

pdf("empirical_results_depth.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(2.5e-6, 2e-4)), bty = "n",
     ylab = expression("Evol. rate"~10^{-5}~"subs/site/year"), 
     yaxt = "n", xaxt = "n", xlab = "")
axis(2, at = log10(c(2.5e-6, 5e-6, 1.5e-5, 5e-5, 1e-4)), 
     labels = c("0.25", "0.5", "1.5", "5.0",  "100"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "Modern\nonly", "Up to 500\nyears BP", "Up to 1,000\nyears BP", "Up to 5,000\nyears BP", ""))

prior <- density(log10(T0_hbv_sample_prior$rate.mean))
posterior <- density(log10(T0_hbv_prior$rate.mean))
polygon(1+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T500_hbv_sample_prior$rate.mean))
posterior <- density(log10(T500_hbv_prior$rate.mean))
polygon(3.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T1000_hbv_sample_prior$rate.mean))
posterior <- density(log10(T1000_hbv_prior$rate.mean))
polygon(6+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T_3020_hbv_sample_prior$rate.mean))
posterior <- density(log10(T_3020_hbv_prior$rate.mean))
polygon(8.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
#
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()



######################################### Plot root height?
#########################################

pdf("empirical_results_tree_height.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
    rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(1, 100000)), bty = "n",
     ylab = "Tree height (years before present)", type = "n",
     xlab = "", xaxt = "n", yaxt = "n")
axis(2, at = c(0, 1, 2, 3, 4, 5), 
     labels = c("0", "10", "100", "1000", "10 Ky", "100 Ky"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "Modern\nonly", "Up to 500\nyears BP", "Up to 1,000\nyears BP", "Up to 5,000\nyears BP", ""))
prior <- density(log10(T0_hbv_sample_prior$tree.height))
posterior <- density(log10(T0_hbv_prior$tree.height))
polygon(1+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/2, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T500_hbv_sample_prior$tree.height))
posterior <- density(log10(T500_hbv_prior$tree.height))
polygon(3.5+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/3.5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T1000_hbv_sample_prior$tree.height))
posterior <- density(log10(T1000_hbv_prior$tree.height))
polygon(6+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/3.6, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T_3020_hbv_sample_prior$tree.height))
posterior <- density(log10(T_3020_hbv_prior$tree.height))
polygon(8.5+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/5, c(posterior$x, posterior$x), col = "blue", border = NA)
dev.off()


#########################################
#########################################



#####
# Now with misleading prior

T0_hbv_prior <- tail(read.table("T0_misleading_priors.log", head = T), 1000)
T0_hbv_sample_prior <- tail(read.table("T0_sample_prior_misleading_priors.log", head = T), 1000)

T500_hbv_prior <- tail(read.table("T500_misleading_priors.log", head = T), 1000)
T500_hbv_sample_prior <- tail(read.table("T500_sample_prior_misleading_priors.log", head = T), 1000)

T1000_hbv_prior <- tail(read.table("T1000_misleading_priors.log", head = T), 1000)
T1000_hbv_sample_prior <- tail(read.table("T1000_sample_prior_misleading_priors.log", head = T), 1000)

T_3020_hbv_prior <- tail(read.table("T_3020_misleading_priors.log", head = T), 1000)
T_3020_hbv_sample_prior <- tail(read.table("T_3020_sample_prior_misleading_priors.log", head = T), 1000)

pdf("empirical_results_depth_misleading_prior.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(5e-7, 2e-4)), bty = "n",
     ylab = expression("Evol. rate"~10^{-5}~"subs/site/year"), 
     yaxt = "n", xaxt = "n", xlab = "")
axis(2, at = log10(c(5e-7, 1.5e-6, 5e-6, 1.5e-5, 5e-5, 1e-4)), 
     labels = c("0.05", "0.15", "0.5", "1.5", "5.0",  "100"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "Modern\nonly", "Up to 500\nyears BP", "Up to 1,000\nyears BP", "Up to 5,000\nyears BP", ""))

prior <- density(log10(T0_hbv_sample_prior$rate.mean))
posterior <- density(log10(T0_hbv_prior$rate.mean))
polygon(1+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T500_hbv_sample_prior$rate.mean))
posterior <- density(log10(T500_hbv_prior$rate.mean))
polygon(3.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T1000_hbv_sample_prior$rate.mean))
posterior <- density(log10(T1000_hbv_prior$rate.mean))
polygon(6+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T_3020_hbv_sample_prior$rate.mean))
posterior <- density(log10(T_3020_hbv_prior$rate.mean))
polygon(8.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
#
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()



############################ Plot root height
T0_hbv_prior <- tail(read.table("T0_misleading_priors.log", head = T), 1000)
T0_hbv_sample_prior <- tail(read.table("T0_sample_prior_misleading_priors.log", head = T), 1000)

T500_hbv_prior <- tail(read.table("T500_misleading_priors.log", head = T), 1000)
T500_hbv_sample_prior <- tail(read.table("T500_sample_prior_misleading_priors.log", head = T), 1000)

T1000_hbv_prior <- tail(read.table("T1000_misleading_priors.log", head = T), 1000)
T1000_hbv_sample_prior <- tail(read.table("T1000_sample_prior_misleading_priors.log", head = T), 1000)

T_3020_hbv_prior <- tail(read.table("T_3020_misleading_priors.log", head = T), 1000)
T_3020_hbv_sample_prior <- tail(read.table("T_3020_sample_prior_misleading_priors.log", head = T), 1000)

pdf("empirical_results_depth_misleading_prior_tree_height.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(phylo_threshold_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(1, 1000000)), bty = "n",
     ylab = "Tree height (years before present)", type = "n",
     xlab = "", xaxt = "n", yaxt = "n")
axis(2, at = c(0, 1, 2, 3, 4, 5, 6, 7), 
     labels = c("0", "10", "100", "1000", "10 Ky", "100 Ky", "1 My", "10 My"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "Modern\nonly", "Up to 500\nyears BP", "Up to 1,000\nyears BP", "Up to 5,000\nyears BP", ""))
prior <- density(log10(T0_hbv_sample_prior$tree.height))
posterior <- density(log10(T0_hbv_prior$tree.height))
polygon(1+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/3.5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T500_hbv_sample_prior$tree.height))
posterior <- density(log10(T500_hbv_prior$tree.height))
polygon(3.5+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/3.5, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T1000_hbv_sample_prior$tree.height))
posterior <- density(log10(T1000_hbv_prior$tree.height))
polygon(6+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/3.6, c(posterior$x, posterior$x), col = "blue", border = NA)

prior <- density(log10(T_3020_hbv_sample_prior$tree.height))
posterior <- density(log10(T_3020_hbv_prior$tree.height))
polygon(8.5+c(prior$y, -prior$y), c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/5, c(posterior$x, posterior$x), col = "blue", border = NA)
#
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()















##############################################
###### Biased sampling results
biased_sampling_hcc_trees <- c("modern95.trees.hcc.tre", "modern50.trees.hcc.tre", 
                               "modern25.trees.hcc.tre", "modern10.trees.hcc.tre")

pdf("empirical_results_biased.pdf", useDingbats = "FALSE", width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
modern95_hbv_prior <- tail(read.table("modern95.log", head = T), 1000)
modern95_hbv_sample_prior <- tail(read.table("modern95_sample_prior.log", head = T), 1000)

modern50_hbv_prior <- tail(read.table("modern50.log", head = T), 1000)
modern50_hbv_sample_prior <- tail(read.table("modern50_sample_prior.log", head = T), 1000)

modern25_hbv_prior <- tail(read.table("modern25.log", head = T), 1000)
modern25_hbv_sample_prior <- tail(read.table("modern25_sample_prior.log", head = T), 1000)

modern10_hbv_prior <- tail(read.table("modern10.log", head = T), 1000)
modern10_hbv_sample_prior <- tail(read.table("modern10_sample_prior.log", head = T), 1000)

par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(2.5e-6, 2e-4)), bty = "n",
     ylab = expression("Evol. rate"~10^{-5}~"subs/site/year"), 
     yaxt = "n", xaxt = "n", xlab = "")
axis(2, at = log10(c(2.5e-6, 5e-6, 1.5e-5, 5e-5, 1e-4)), 
     labels = c("0.25", "0.5", "1.5", "5.0",  "100"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "95% modern", "50% modern", "25% modern", "10% modern", ""))

prior <- density(log10(modern95_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern95_hbv_prior$rate.mean))
polygon(1+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern50_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern50_hbv_prior$rate.mean))
polygon(3.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern25_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern25_hbv_prior$rate.mean))
polygon(6+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern10_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern10_hbv_prior$rate.mean))
polygon(8.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()

############################ Biased sampling results - Root height
biased_sampling_hcc_trees <- c("modern95.trees.hcc.tre", "modern50.trees.hcc.tre", 
                               "modern25.trees.hcc.tre", "modern10.trees.hcc.tre")

pdf("empirical_results_biased_root_height.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
modern95_hbv_prior <- tail(read.table("modern95.log", head = T), 1000)
modern95_hbv_sample_prior <- tail(read.table("modern95_sample_prior.log", head = T), 1000)

modern50_hbv_prior <- tail(read.table("modern50.log", head = T), 1000)
modern50_hbv_sample_prior <- tail(read.table("modern50_sample_prior.log", head = T), 1000)

modern25_hbv_prior <- tail(read.table("modern25.log", head = T), 1000)
modern25_hbv_sample_prior <- tail(read.table("modern25_sample_prior.log", head = T), 1000)

modern10_hbv_prior <- tail(read.table("modern10.log", head = T), 1000)
modern10_hbv_sample_prior <- tail(read.table("modern10_sample_prior.log", head = T), 1000)

par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(1, 100000)), bty = "n",
     ylab = "Tree height (years before present)", type = "n",
     xlab = "", xaxt = "n", yaxt = "n")
axis(2, at = c(0, 1, 2, 3, 4, 5), 
     labels = c("0", "10", "100", "1000", "10 Ky", "100 Ky"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "95% modern", "50% modern", "25% modern", "10% modern", ""))

prior <- density(log10(modern95_hbv_sample_prior$tree.height))
posterior <- density(log10(modern95_hbv_prior$tree.height))
polygon(1+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern50_hbv_sample_prior$tree.height))
posterior <- density(log10(modern50_hbv_prior$tree.height))
polygon(3.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern25_hbv_sample_prior$tree.height))
posterior <- density(log10(modern25_hbv_prior$tree.height))
polygon(6+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern10_hbv_sample_prior$tree.height))
posterior <- density(log10(modern10_hbv_prior$tree.height))
polygon(8.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()






##############################################
###### Biased sampling results- Misleading prior
biased_sampling_hcc_trees <- c("modern95.trees.hcc.tre", "modern50.trees.hcc.tre", 
                               "modern25.trees.hcc.tre", "modern10.trees.hcc.tre")

pdf("empirical_results_biased_misleading_prior.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
modern95_hbv_prior <- tail(read.table("modern95_misleading_priors.log", head = T), 1000)
modern95_hbv_sample_prior <- tail(read.table("modern95_sample_prior_misleading_priors.log", head = T), 1000)

modern50_hbv_prior <- tail(read.table("modern50_misleading_priors.log", head = T), 1000)
modern50_hbv_sample_prior <- tail(read.table("modern50_sample_prior_misleading_priors.log", head = T), 1000)

modern25_hbv_prior <- tail(read.table("modern25_misleading_priors.log", head = T), 1000)
modern25_hbv_sample_prior <- tail(read.table("modern25_sample_prior_misleading_priors.log", head = T), 1000)

modern10_hbv_prior <- tail(read.table("modern10_misleading_priors.log", head = T), 1000)
modern10_hbv_sample_prior <- tail(read.table("modern10_sample_prior_misleading_priors.log", head = T), 1000)

par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(5e-7, 2e-4)), bty = "n",
     ylab = expression("Evol. rate"~10^{-5}~"subs/site/year"), 
     yaxt = "n", xaxt = "n", xlab = "")
axis(2, at = log10(c(5e-7, 1.5e-6, 5e-6, 1.5e-5, 5e-5, 1e-4)), 
     labels = c("0.05", "0.15", "0.5", "1.5", "5.0",  "100"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "95% modern", "50% modern", "25% modern", "10% modern", ""))

prior <- density(log10(modern95_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern95_hbv_prior$rate.mean))
polygon(1+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern50_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern50_hbv_prior$rate.mean))
polygon(3.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern25_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern25_hbv_prior$rate.mean))
polygon(6+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern10_hbv_sample_prior$rate.mean))
posterior <- density(log10(modern10_hbv_prior$rate.mean))
polygon(8.5+c(prior$y, -prior$y)/5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()




##############################################
###### Biased sampling results- Misleading prior -root height
biased_sampling_hcc_trees <- c("modern95.trees.hcc.tre", "modern50.trees.hcc.tre", 
                               "modern25.trees.hcc.tre", "modern10.trees.hcc.tre")

pdf("empirical_results_biased_misleading_prior_root_height.pdf", useDingbats = FALSE, width = 8, height = 9)
par(mfrow = c(2, 1), mar = c(0, 4.4, 1, 1))
plot(1, 1, type = "n", xlim = c(0, 500), ylim = log10(c(1, 16000)), bty = "n", 
     xaxt = "n", xlab = "", ylab = "Time (years BP)", yaxt = "n")
axis(2, at = log10(c(1, 20, 1000, 10000, 20000)), 
     labels = c("0", "20", "1K", "10K", "20K"), las = 2)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[1])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[2])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 120)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[3])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 240)
tr_temp <- read.tree(text = write.tree(read.nexus(biased_sampling_hcc_trees[4])))
plot.tree.lines(tr_temp,  plot.new = F, line.type = "l", 
                rotation.angle = pi*3/2, log.scale = T,    x.offset = 360)
#
modern95_hbv_prior <- tail(read.table("modern95_misleading_priors.log", head = T), 1000)
modern95_hbv_sample_prior <- tail(read.table("modern95_sample_prior_misleading_priors.log", head = T), 1000)

modern50_hbv_prior <- tail(read.table("modern50_misleading_priors.log", head = T), 1000)
modern50_hbv_sample_prior <- tail(read.table("modern50_sample_prior_misleading_priors.log", head = T), 1000)

modern25_hbv_prior <- tail(read.table("modern25_misleading_priors.log", head = T), 1000)
modern25_hbv_sample_prior <- tail(read.table("modern25_sample_prior_misleading_priors.log", head = T), 1000)

modern10_hbv_prior <- tail(read.table("modern10_misleading_priors.log", head = T), 1000)
modern10_hbv_sample_prior <- tail(read.table("modern10_sample_prior_misleading_priors.log", head = T), 1000)

par(mar = c(4, 4.4, 1, 1))
plot(1, 1, xlim = c(0, 10), ylim = log10(c(1, 1000000)), bty = "n",
     ylab = "Tree height (years before present)", type = "n",
     xlab = "", xaxt = "n", yaxt = "n")
axis(2, at = c(0, 1, 2, 3, 4, 5, 6, 7), 
     labels = c("0", "10", "100", "1000", "10 Ky", "100 Ky", "1 My", "10 My"), las = 2)
axis(1, at = c(0, 1, 3.5, 6, 8.5, 10), padj = 0.7,
     labels = c("", "95% modern", "50% modern", "25% modern", "10% modern", ""))

prior <- density(log10(modern95_hbv_sample_prior$tree.height))
posterior <- density(log10(modern95_hbv_prior$tree.height))
polygon(1+c(prior$y, -prior$y)/2.5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(1+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern50_hbv_sample_prior$tree.height))
posterior <- density(log10(modern50_hbv_prior$tree.height))
polygon(3.5+c(prior$y, -prior$y)/2.5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(3.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern25_hbv_sample_prior$tree.height))
posterior <- density(log10(modern25_hbv_prior$tree.height))
polygon(6+c(prior$y, -prior$y)/2.5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(6+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
prior <- density(log10(modern10_hbv_sample_prior$tree.height))
posterior <- density(log10(modern10_hbv_prior$tree.height))
polygon(8.5+c(prior$y, -prior$y)/2.5, c(prior$x, prior$x), col = "orange", border = NA)
polygon(8.5+c(posterior$y, -posterior$y)/10, c(posterior$x, posterior$x), col = "blue", border = NA)
lines(c(-1, 10), log10(c(1.5e-5, 1.5e-5)), lty = 2)
dev.off()

