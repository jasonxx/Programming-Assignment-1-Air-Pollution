library(fields);  # Neeeds this for image.plot
library(stringr); # Need for extracting residue names using regular expression

old.par <- par(no.readonly=T);

# This script is meant to work under Linux. 
# Windows may use a different convention for the "\" in directory paths.

# The script assumes the following directories in the working directory to exist:
# out : plots and data will be saved here
# list : contains the Sparky list files in this format (without the initial "#"):
#Assignment         w1         w2        S/N    lw1 (hz)   lw2 (hz) 
#
#A3N-HN    125.706      8.427        324       34.6       26.7 
#A4N-HN    123.372      8.031        584       30.7       23.1 
#A5N-HN    123.628      8.060        581       33.6       25.0 
# ...
# Note that the script will take care or removing the initial two lines of the list files.


# The script will read in all files in the "list" folder that have a ".list" suffix, like "filename.list".
# Make sure they all belong to the same project!

# MAKE ALL CHANGES HERE
# ---------------------------------------------------------------------------------------
# Set your working directory. 
setwd('C:/Users/USER/Documents/concise');
methyl    <- 0;      # Backbone amide (0) or side chain methyl group (1)?
refere    <- "pki";  # "0" is spectra were already referenced in Sparky, "filename" of the reference if not.
l.state   <- "apo";  # "Left" state. These two have to match the filename.list
r.state   <- "pki";  # "Right" state.
pc.states <- 0; # States used to identify the linearity. If "0" all states are used. Otherwise specify
                # the states, eg: c("amppnp", "amppnp_pki", "amppnp_pln", "apo")
#pc.states <- c("amppnp", "amppnp_pki", "amppnp_pln", "amppnp_r14del", "apo");
ch.cutoff <- 0.97; # The CHESCA cutoff, see PNAS Melacini. For sample size 5-6 0.97 should be OK
                   # For larger sample size, maybe it could be reduced.
                   # Increase to 0.99 for sample of 4 
# ---------------------------------------------------------------------------------------
# END OF ALL CHANGES

# CONSTANTS
# ---------------------------------------------------------------------------------------
lin.fil    <- 3.0; # SD(PC1)/SD(PC2) linearity filter threshold;
ppm.fil    <- 0.05 - methyl*0.02; # Min ppm distance below which discard. It's 0.05 for HSQC and 0.03 for HMQC
regex      <- "([[:alpha:]]+)([[:digit:]]+)([[:alnum:]]+)-([[:alnum:]]+)";
s.scale    <- 0.146 * methyl + 0.154; # Heavy atom scaling factor. S=N 0.154;   S=C 0.300
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
# ---------------------------------------------------------------------------------------

# FUNCTIONS
# ---------------------------------------------------------------------------------------
# SPEC_DIST
# -----------------------------
# Calculates the distance in ppm between two spectra
spec_dist <- function(par) {
  # B is the reference
	dn <- par[1];
	dh <- par[2];
	am <- a;
	bm <- b;
	am[, 1] <- am[, 1] + dn;
	am[, 2] <- am[, 2] + dh;
	d <- sum(abs(am - b));
	return(d);
}
# -----------------------------

# FIND LINEARITY THROUGH PCA 
# -----------------------------------------------------------------------------
linearity <- function(dxy, zero=1, one=1, ref=1) {
  # Reads in X and Y, does PCA and reorders putting Left before Right
  # It then scales the scores to normalize to one standard deviation
  # Returns, all in one line, the following data:
  # 1) xtreme: "0" if Left and Right are not at the extremes, "1" otherwise
  # 2:3) sd.pc.1 and sd.pc.2
  # 4) x.pc.range: the maximum minus the minimum along x, unscaled
  # 5:5+n.row-1) pc1.scores, unscaled
  # 5+n.row:5+2*n.row-1) pc.2.scores, unscaled
  # 5+2*n.row:5+3*n.row-1) pc.1.scores, scaled
  n.row = length(dxy)/2;
  xtra <- (1:n.row)[!(1:n.row %in% ref)];
  xy.all <- matrix(dxy, nrow=n.row, ncol=2, byrow=T);
  xy <- xy.all[ref, ];
  xy.xtra <- matrix(xy.all[xtra, ], nrow=length(xtra), ncol=2);
  pc <- princomp(xy);
  xy.pc <- pc$scores;
  xtra.pc <- predict(pc, xy.xtra);
  xy.pc <- rbind(xy.pc, xtra.pc);
  flip <- (xy.pc[one, 1] - xy.pc[zero, 1])/abs(xy.pc[one, 1] - xy.pc[zero, 1]);
  x.pc <- flip*xy.pc[ , 1];
  y.pc <- flip*xy.pc[ , 2];
  x.pc.range <- abs(x.pc[zero] - x.pc[one]);
  x.sc <- x.pc / pc$sdev[1];
  xtreme <- ifelse((x.pc.range < (max(x.pc) - min(x.pc))) ,0 ,1);
  linear.summ <- as.numeric(c(xtreme, pc$sdev, x.pc.range, x.pc, y.pc, x.sc));
}
# -----------------------------------------------------------------------------

# READ THE DATA
# ---------------------------------------------------------------------------------------
list.states  <- vector();
df           <- data.frame();
list.files   <- (Sys.glob("list/*.list"));
for (file in list.files) {
  state         <- substr(file, 6, nchar(file)-5);  # removes list/ and .list
  list.states   <- append(list.states, state); 
  temp          <- read.table(file, header=F, skip=2, blank.lines.skip=T, as.is=T);
  temp.df       <- data.frame(cbind(state, str_match(temp$V1, regex), temp$V2, temp$V3), stringsAsFactors=F);
  df            <- rbind(df, temp.df);
}
# ---------------------------------------------------------------------------------------

# CLEAN THE DATA
# ---------------------------------------------------------------------------------------
n.states  <- length(list.states);
colnames(df)  <- c("state", "label", "resn", "resi", "Sname", "Iname", "Sppm", "Ippm");
df$state  <- as.factor(df$state);
df$label  <- as.factor(df$label);
df$resi   <- as.numeric(df$resi);
df$Sppm   <- as.numeric(df$Sppm);
df$Ippm   <- as.numeric(df$Ippm);
# Only keep the "labels" that are common to all states
suppressWarnings(temp.df   <- Reduce(function(...) merge(..., all=F, by="label"), split(df, df$state)));
# Not the most elegant way, but just keep some basic info and the S and I ppm.
keep      <- sort(c(1, 3:6, 1:n.states * 7, 1:n.states * 7 + 1));
df.comm   <- temp.df[, keep];       # In df.comm are kept only the labels common to all states
st.names  <- temp.df[1, (1:n.states-1)*7 + 2]; # Should be the same as list.states, but just in case they got scrambled...
col.names <- c("label", "resn", "resi", "Sname", "Iname", sapply(st.names, function(x) c(paste(x, ".Sppm", sep=""), paste(x, ".Ippm", sep=""))));
colnames(df.comm) <- col.names;
df.comm   <- df.comm[order(df.comm$resi, df.comm$label), ]; # Order first by resid and then by label
s.col     <- (1:n.states-1)*2 + 6; # S nuclues columns
df.scal   <- df.comm;
df.scal[, s.col] <- df.scal[, s.col] * s.scale;       # In df.scal we have scaled the heavy nucleus ppm
s.mat     <- data.matrix(df.scal[, 6:ncol(df.scal)]); # In s.mat we keep only the chem shift
s.orig    <- s.mat;
si        <- seq(1, ncol(s.mat), 2);
resi      <- df.comm$resi;
resn      <- df.comm$resn;
label     <- df.comm$label;
Sname     <- df.comm$Sname;
# ---------------------------------------------------------------------------------------

# RECENTER SPECTRA
# ---------------------------------------------------------------------------------------
if (refere != 0) {
  ref <- which(list.states == refere, arr.ind=T);  
  b <- s.mat[ , c(si[ref], si[ref]+1)];
  for (i in si) {
   	a <- s.mat[ , i:(i+1)];
	  opt <- optim(c(0, 0), spec_dist, control=list(abstol=0.0001));
	  dnh <- opt$par;
	  if (abs(dnh[1]) < 0.001) {dnh[1] <- 0}
	  if (abs(dnh[2]) < 0.001) {dnh[2] <- 0}
	  s.mat[ , i] <- s.mat[ , i] + dnh[1];
	  s.mat[ , (i+1)] <- s.mat[ , (i+1)] + dnh[2];
  }
}
# ---------------------------------------------------------------------------------------

# LINEARITY
# -----------------------------------------------------------------------------
n.left   <- which (list.states %in% l.state, arr.ind=T);
n.right  <- which (list.states %in% r.state, arr.ind=T);
s   <- 1:n.states;
si  <- (s - 1)*2 + 1;
if (length(pc.states) == 1) {
  if (pc.states == 0) {
    pc.states <- list.states;
  }
}
# -----------------------------------------------------------------------------

# DATA ANALYSIS AND PLOT
# -----------------------------------------------------------------------------
s.lincol <- 1:n.states+4+n.states*2;
s.ref    <- which(list.states %in% pc.states, arr.ind=T);
m <- t(apply(s.mat, 1, linearity, zero=n.left, one=n.right, ref=s.ref));

# No Filter
sc.m <- colMeans(m[, s.lincol]);
sc.sd <- apply(m[, s.lincol], 2, sd);
xx   <- seq(-3.0, 3.0, 0.001);
ymax <- 0;
for (i in 1:n.states) {ymax <- max(c(ymax, dnorm(xx, mean=sc.m[i], sd=sc.sd[i])))}
ymax <- round(ymax + 0.1, 1);
pdf('out/concise_nofilter.pdf');
par(cex=1.5);
plot (-999, 999, xlim=range(xx), ylim=c(0, ymax), xlab=paste(l.state, r.state, sep=" -> "), ylab='Probability Density', cex=2.0);
for (i in 1:n.states) {
  lines(xx, dnorm(xx, mean=sc.m[i], sd=sc.sd[i]), t='l', col=jet.colors(n.states)[i], lwd=3);
  text(-3, ymax-i*0.05, list.states[i], col=jet.colors(n.states)[i], pos=4, cex=0.6);
}
dev.off();

# Filter
reduced.set   <- which(m[,2]/m[,3] > lin.fil & m[,4] > ppm.fil, arr.ind=T);
sc.m.reduced  <- colMeans(m[reduced.set, s.lincol]);
sc.sd.reduced <- apply(m[reduced.set, s.lincol], 2, sd);
ymax <- 0;
for (i in 1:n.states) {ymax <- max(c(ymax, dnorm(xx, mean=sc.m.reduced[i], sd=sc.sd.reduced[i])))}
ymax <- round(ymax + 0.1, 1);
pdf('out/concise_filter.pdf');
par(cex=1.5);
plot (-999, 999, xlim=range(xx), ylim=c(0, ymax), xlab=paste(l.state, r.state, sep=" -> "), ylab='Probability Density', cex=2.0);
for (i in 1:n.states) {
  lines(xx, dnorm(xx, mean=sc.m.reduced[i], sd=sc.sd.reduced[i]), t='l', col=jet.colors(n.states)[i], lwd=3);
  text(-3, ymax-i*0.05, list.states[i], col=jet.colors(n.states)[i], pos=4, cex=0.6);
}
dev.off();

# Results
cat("NO-Filter\n# ------\nResidue list\n", file="out/residue_list.dat");
cat(resi, file="out/residue_list.dat", sep="+", append=T);
cat("\n", file="out/residue_list.dat", append=T);
cat(as.character(label), file="out/residue_list.dat", sep=" ", append=T);
cat("\n", file="out/residue_list.dat", append=T);
cat("# ------\n\nFilter\n# ------\nResidue list\n", file="out/residue_list.dat", append=T);
cat(resi[reduced.set], file="out/residue_list.dat", sep="+", append=T);
cat("\n", file="out/residue_list.dat", append=T);
cat(as.character(label[reduced.set]), file="out/residue_list.dat", sep=" ", append=T);
cat("\n# -----\n\n", file="out/residue_list.dat", append=T);

sc.m.ord <- sort(sc.m, index.return=T)$ix;
sc.m.red.ord <- sort(sc.m.reduced, index.return=T)$ix;
tab <- data.frame(rbind(
       list.states[sc.m.ord], 
       format(sc.m[sc.m.ord], digits=2, nsmall=2), 
       format(sc.sd[sc.m.ord], digits=2, nsmall=2),
       list.states[sc.m.red.ord], 
       format(sc.m.reduced[sc.m.red.ord], digits=2, nsmall=2), 
       format(sc.sd.reduced[sc.m.red.ord], digits=2, nsmall=2)),
       row.names=c("NoFilt State", "NoFilt Mean", "NoFilt SD", "Filt State", "Filt Mean", "Filt SD"));
capture.output(print(tab, print.gap=3), file="out/concise.dat");

# CHESCA (Performed only on REDUCED SET)
# ---------------------------------------
colors <- jet.colors(5);
h      <- 1 - ch.cutoff;
l      <- length(resi);
reds   <- resi[reduced.set];
t.tag.char <- label[reduced.set];
len <- length(t.tag.char);
M <- m[reduced.set, s.lincol];
R <- cor(t(M));
Rabs <- abs(R);
img.len <- max(resi); # Only used for NH
res     <- 1:img.len;    # Only used for NH
if (methyl) {
  image.plot(1:len, 1:len, R, col=colors, zlim=c(0.9, 1), nlevel=length(colors), xlab='PKA ResID', ylab='PKA ResID', axes=T);
  pdf('out/chesca_correlations.pdf');
    par(mar=c(6.1,6.1,5.1,5.1));
    image(1:len, 1:len, abs(R), axes=F,col='transparent', xlab="", ylab="");
    par(cex.axis=0.75, las=2);
    axis(2, at=1:len, labels=t.tag.char);
    axis(1, at=1:len, labels=t.tag.char);
    par(cex.axis=1.0, las=2);
    image.plot(1:len, 1:len, abs(R), col=colors, nlevel=length(colors), zlim=c(0.9, 1), axes=F, legend.mar=4.1, add=T);
    rect(0, 0, len+0.5, len+0.5);
  dev.off();
} else {
  R.img <- array(0.0, dim=c(img.len, img.len));
  R.img[reds, reds] <- Rabs;
  pdf('out/chesca_correlations.pdf');
    par(mar=c(5, 4, 5.2, 2), lab=c(10, 10, 10));
    image.plot(res, res, R.img, col=colors, zlim=c(0.9, 1), nlevel=length(colors), xlab='PKA ResID', ylab='PKA ResID', axes=T);
  dev.off();
}

if (methyl) {
  rac.names <- t.tag.char;
} else {
  rac.names <- reds;
}

diss <- as.dist(1 - Rabs);
RAC <- hclust(diss, method='single');
pdf(file='out/chesca_dendrogram.pdf', height=7, width=9);
  plot(RAC, labels=rac.names, xlab="ResID");
  rac.plot <- rect.hclust(RAC, h=h, border=2);
dev.off();
sum.rac <- summary(rac.plot);
rac.len <- as.numeric(sum.rac[,1]);
rac.1.2 <- sort(rac.len, decreasing=T)[c(1,2)];
rac.1   <- which(rac.len %in% rac.1.2[1], arr.ind=T);
rac.2   <- which(rac.len %in% rac.1.2[2], arr.ind=T);
cl.1.inv <- as.numeric(rac.plot[[rac.1]]);
cl.1 <- rac.names[cl.1.inv]; 

Rcut <- Rabs;
Rcut[Rcut < ch.cutoff] <- NA;
pairs <- which(Rabs >= Rcut, arr.ind=T);
pairs.resi.1 <- pairs[,1];
pairs.resi.2 <- pairs[,2];
if (methyl) {
  min.1 <- min(cl.1.inv);
  max.1 <- max(cl.1.inv);
} else {
  min.1 <- min(cl.1);
  max.1 <- max(cl.1);
}

pdf(file='out/chesca_matrix.pdf', height=7, width=7);
  par(cex.axis=0.75, las=2, cex=1.0, mar=c(6.0,6.0,4.0,4.0));
  if (methyl) {
    plot(pairs.resi.1, pairs.resi.2, pch=19, xlab='', ylab='', axes=F);
    axis(2, at=1:len, labels=t.tag.char);
    axis(1, at=1:len, labels=t.tag.char);
    box();
    for (res in cl.1.inv) {
      segments(x0=min.1, y0=res, x1=max.1, y1=res, col='blue', lwd=1.5);
      segments(y0=min.1, x0=res, y1=max.1, x1=res, col='blue', lwd=1.5);
    }
  } else {
    plot(reds[pairs.resi.1], reds[pairs.resi.2], pch=19, xlab='PKA ResID', ylab='PKA ResID', xlim=c(0, img.len), ylim=c(0, img.len));
    for (ires in cl.1) {
      segments(x0=min.1, y0=ires, x1=max.1, y1=ires, col='blue', lwd=1.5);
      segments(y0=min.1, x0=ires, y1=max.1, x1=ires, col='blue', lwd=1.5);
    }
  }
dev.off();

if (methyl) {
  labels.df <- data.frame(str_match(cl.1, regex), stringsAsFactors=F);
  pymol.tmp <- apply( labels.df[ , c(3, 4)] , 1 , paste , collapse = " & n. ");
  pymol.sel <- paste("sel cluster, i. ", paste(pymol.tmp, collapse=" + i. "), sep="");
  
  cat("CHESCA Main Cluster\n# ------\nUnique Residue list\n", file="out/chesca_list.dat");
  cat(as.character(unique(labels.df[, 3])), file="out/chesca_list.dat", sep=" ", append=T);
  cat("\nSpecific atoms list\n", file="out/chesca_list.dat", sep=" ", append=T);
  cat(as.character(cl.1), file="out/chesca_list.dat", sep=" ", append=T);
  cat("\n", file="out/residue_list.dat", append=T);
  cat("\nPyMOL Selection\n", file="out/chesca_list.dat", sep=" ", append=T);
  cat(as.character(pymol.sel), file="out/chesca_list.dat", sep=" ", append=T);
} else {
  pymol.sel <- paste("sel cluster, n. CA & i. ", paste(cl.1, collapse="+"), sep="");
  
  cat("CHESCA Main Cluster\n# ------\nUnique Residue list\n", file="out/chesca_list.dat");
  cat(as.character(cl.1), file="out/chesca_list.dat", sep=" ", append=T);
  cat("\n", file="out/residue_list.dat", append=T);
  cat("\nPyMOL Selection\n", file="out/chesca_list.dat", sep=" ", append=T);
  cat(as.character(pymol.sel), file="out/chesca_list.dat", sep=" ", append=T);
}


