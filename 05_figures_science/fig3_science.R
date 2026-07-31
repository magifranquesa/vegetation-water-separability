#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_figures/fig3.R
#
# Project: Aridity and timescale bound the separability of vegetation water stress
#
# Figure 3: vegetation-hydroclimate correlation strength and the composition of
# hydroclimatic control along the aridity gradient, on a continuous AI axis.
#   Top (a):    area-weighted mean |rho*| in fine AI bins (0.02), loess-smoothed
#               (span 0.4), per indicator -- binned means plus a local smoother
#               (not a GAM; no threshold is read off this curve). All vegetated
#               cells regardless of significance (strength; extent is Fig. 4).
#               One row = panel "a"; the four columns are month strips.
#   Bottom (b): area-weighted share of the significant domain under each control
#               vs AI (fine bins, loess), with the soil->joint crossover marked.
#               One row = panel "b"; the four columns are month strips.
#
# Inputs (READ-ONLY):
#   outputs/intermediate/absmax_spearman/abs_max_correlation_kndvi_{Ep,Et,ED,SMrz,SMs}.nc
#     (variable abs_max_correlation; panel a)
#   outputs/varpart_global/varpart_signif_global_2blocks.nc   (total_r2; panel b)
#   outputs/varpart_global/fdr_adjusted_pvalues.nc            (p_full/soil/demand_adj; panel b)
#   data/processed/aridity/ai_1982_2022_period.nc             (aridity index)
#   outputs/intermediate/vegetation_mask_c1.tif
#
# Output:
#   outputs/figures_science/fig3_plateau_transition_science.tif
#   (PDF output commented out)
# ==============================================================================

suppressPackageStartupMessages({ library(terra); library(ggplot2); library(patchwork); library(scales) })
config_file <- file.path("R","config.R"); if (file.exists(config_file)) source(config_file)

dir_absmax <- file.path("outputs","intermediate","absmax_spearman")
file_nc  <- file.path("outputs","varpart_global","varpart_signif_global_2blocks.nc")
file_adj <- file.path("outputs","varpart_global","fdr_adjusted_pvalues.nc")
file_ai  <- if (exists("paths")&&!is.null(paths$aridity_index)) paths$aridity_index else file.path("data","processed","aridity","ai_1982_2022_period.nc")
file_veg <- if (exists("paths")&&!is.null(paths$veg_mask_c1)) paths$veg_mask_c1 else file.path("outputs","intermediate","vegetation_mask_c1.tif")
figures_dir <- file.path("outputs","figures_science")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
out_tif <- file.path(figures_dir,"fig3_plateau_transition_science.tif")
# out_pdf <- sub("\\.tif$", ".pdf", out_tif)   # PDF output disabled

indicators <- c("Ep","Et","ED","SMrz","SMs"); ind_lab <- c(Ep="AED",Et="Et",ED="ED",SMrz="SMrz",SMs="SMs")
ind_colours <- c(AED="#D98C00",Et="#4DAF4A",ED="#9E1F63",SMrz="#084C8D",SMs="#6BAED6")
ctrl_levels <- c("Soil moisture","Joint (non-separable)","Atmospheric demand")
ctrl_colours <- c("Soil moisture"="#2166AC","Joint (non-separable)"="#5E6B7A","Atmospheric demand"="#D98C00")
months_index <- c(1,4,7,10); months_names <- c("January","April","July","October")
ai_thresholds <- c(0.03,0.20,0.50,0.65); FV <- -9999; ALPHA <- 0.05; XLIM <- c(0,0.9); SPAN <- 0.4
font_family <- "sans"; base_size <- 7; lab_size <- 8; fig_width_mm <- 183; fig_height_mm <- 92; fig_dpi <- 300
set_na <- function(v){ v[abs(v-FV)<1] <- NA; v }
bw <- 0.02; brks <- seq(0,1.6,by=bw); mids <- head(brks,-1)+bw/2

# ---- TOP: area-weighted mean |rho*| in fine AI bins ----
message("Top: area-weighted mean |rho*| in fine AI bins ...")
ref <- rast(file.path(dir_absmax,"abs_max_correlation_kndvi_Ep.nc"), subds="abs_max_correlation")[[1]]
if (is.na(crs(ref))||crs(ref)=="") crs(ref) <- "EPSG:4326"
area_v <- as.vector(values(cellSize(ref,unit="km",mask=FALSE)))
veg <- as.vector(values(resample(rast(file_veg),ref,method="near"))); veg <- !is.na(veg)&veg>0
ai <- rast(file_ai)[[1]]; if (is.na(crs(ai))) crs(ai) <- "EPSG:4326"; if (xmax(ai)>180+1e-6) ai <- rotate(ai)
ai_v <- as.vector(values(resample(ai,ref,method="bilinear"))); bin_v <- cut(ai_v,breaks=brks,right=FALSE,include.lowest=TRUE)
top_rows <- list()
for (v in indicators){ rc <- rast(file.path(dir_absmax,paste0("abs_max_correlation_kndvi_",v,".nc")), subds="abs_max_correlation")
  for (mi in seq_along(months_index)){ rho <- abs(as.vector(values(rc[[months_index[mi]]])))
    s <- veg & is.finite(rho) & !is.na(bin_v)
    num <- tapply(rho[s]*area_v[s], bin_v[s], sum); den <- tapply(area_v[s], bin_v[s], sum)
    mb <- as.numeric(num)/as.numeric(den); ab <- as.numeric(den)
    ok <- is.finite(mb) & mids>=XLIM[1] & mids<=XLIM[2] & is.finite(ab) & ab>0
    top_rows[[length(top_rows)+1L]] <- data.frame(indicator=ind_lab[[v]], month=months_names[mi], AI=mids[ok], mean_abs_rho=mb[ok], w=ab[ok]) } }
top <- do.call(rbind,top_rows); top$indicator <- factor(top$indicator,levels=c("AED","Et","ED","SMrz","SMs")); top$month <- factor(top$month,levels=months_names)

# ---- BOTTOM: control share vs AI (fine bins) ----
message("Bottom: control share vs AI ...")
r_tt <- rast(file_nc,subds="total_r2"); r_pf <- rast(file_adj,subds="p_full_adj"); r_ps <- rast(file_adj,subds="p_soil_adj"); r_pd <- rast(file_adj,subds="p_demand_adj")
if (is.na(crs(r_tt))||crs(r_tt)=="") crs(r_tt) <- "EPSG:4326"
area2 <- as.vector(values(cellSize(r_tt[[1]],unit="km",mask=FALSE))); ai2 <- as.vector(values(resample(ai,r_tt[[1]],method="bilinear")))
classify_m <- function(m){ tt<-set_na(values(r_tt[[m]])[,1]); pf<-set_na(values(r_pf[[m]])[,1]); ps<-set_na(values(r_ps[[m]])[,1]); pd<-set_na(values(r_pd[[m]])[,1])
  valid<-!is.na(tt); fs<-valid&!is.na(pf)&pf<ALPHA; ss<-!is.na(ps)&ps<ALPHA; ds<-!is.na(pd)&pd<ALPHA
  code<-rep(NA_integer_,length(tt)); code[fs&ss&!ds]<-2L; code[fs&!ss&ds]<-3L; code[fs&!ss&!ds]<-4L; code[fs&ss&ds]<-5L; code }
comp_rows<-list(); cross<-setNames(rep(NA,4),months_names)
for (mi in seq_along(months_index)){ code<-classify_m(months_index[mi]); keep<-!is.na(code)&is.finite(ai2)
  cc<-code[keep]; aa<-area2[keep]; xi<-ai2[keep]; bin<-cut(xi,breaks=brks,right=FALSE,include.lowest=TRUE)
  asum<-function(k){ s<-tapply(aa[cc==k],bin[cc==k],sum); s[is.na(s)]<-0; as.numeric(s) }
  m2<-cbind(soil=asum(2L),joint=asum(4L),dem=asum(3L)); den<-rowSums(m2); ok<-den>0&mids>=XLIM[1]&mids<=XLIM[2]
  comp_rows[[length(comp_rows)+1L]]<-data.frame(month=months_names[mi], AI=rep(mids[ok],3),
    frac=100*c(m2[ok,"soil"],m2[ok,"joint"],m2[ok,"dem"])/rep(den[ok],3),
    control=factor(rep(ctrl_levels,each=sum(ok)),levels=ctrl_levels))
  sf<-100*m2[,"soil"]/ifelse(den>0,den,NA); jf<-100*m2[,"joint"]/ifelse(den>0,den,NA); d<-jf-sf
  cx<-which(d>0&c(FALSE,head(d,-1)<=0)&den>0); cross[mi]<-if(length(cx)) mids[cx[1]] else NA }
comp<-do.call(rbind,comp_rows); comp$month<-factor(comp$month,levels=months_names)

# ---- panels ----
# tags OUTSIDE the top-left corner + top margin so a-h never clip nor hit '100';
# legends built manually (below) so the 8 panels sit in ONE flat layout and their
# columns align top-to-bottom (nested wrap_plots+collect broke the alignment).
theme_p <- function() theme_minimal(base_size=base_size,base_family=font_family) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=base_size,margin=margin(b=1,unit="mm")),
        plot.tag=element_text(face="bold",size=lab_size), plot.tag.position=c(0.01,1.10),
        panel.grid.major=element_line(linewidth=0.25,colour="#E8E8E8"), panel.grid.minor=element_blank(),
        axis.ticks=element_line(linewidth=0.25,colour="#333333"), axis.ticks.length=unit(1,"mm"),
        axis.title=element_text(size=base_size,colour="black"), axis.text=element_text(size=base_size-0.5,colour="#3A3A3A"),
        legend.position="none", plot.margin=margin(3.2,3,1,3,unit="mm"))
# row-level panel tags: only the first column carries a letter (month is the strip
# title), so the top row is a single panel "a" and the bottom row a single panel "b"
ptop<-function(mi){ mn<-months_names[mi]; d<-top[top$month==mn,]
  ggplot(d, aes(AI, mean_abs_rho, colour=indicator, weight=w)) +
    geom_vline(xintercept=ai_thresholds, linetype="dashed", linewidth=0.30, colour="grey60") +
    geom_smooth(method="loess", span=SPAN, se=FALSE, linewidth=0.4) +
    scale_colour_manual(values=ind_colours,name=NULL) +
    scale_x_continuous(breaks=seq(0,0.8,0.2), labels=number_format(accuracy=0.1)) +
    scale_y_continuous(breaks=c(0.2,0.3,0.4,0.5), labels=number_format(accuracy=0.1)) +
    coord_cartesian(xlim=XLIM, ylim=c(0.2,0.55)) +
    labs(title=mn, tag=if(mi==1) "A" else NULL, y=if(mi==1) expression("mean |"*rho*"*|") else NULL) +
    theme_p() + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) }
pbot<-function(mi){ mn<-months_names[mi]; d<-comp[comp$month==mn,]; cx<-cross[mn]
  p<-ggplot(d, aes(AI, frac, colour=control)) + geom_vline(xintercept=ai_thresholds, linetype="dashed", linewidth=0.30, colour="grey60")
  if(!is.na(cx)) p<-p+geom_vline(xintercept=cx, linetype="dotted", linewidth=0.30, colour="#1A1A1A")
  p<-p+geom_smooth(method="loess", span=SPAN, se=FALSE, linewidth=0.4) + scale_colour_manual(values=ctrl_colours,name=NULL) +
    scale_x_continuous(breaks=seq(0,0.8,0.2), labels=number_format(accuracy=0.1)) + scale_y_continuous(breaks=c(0,25,50,75,100)) +
    coord_cartesian(xlim=XLIM, ylim=c(0,100)) + labs(tag=if(mi==1) "B" else NULL, x=NULL, y=if(mi==1) "Significant area (%)" else NULL) +
    theme_p()
  if(!is.na(cx)) p<-p+annotate("text",x=cx,y=99,hjust=-0.1,vjust=1,label=sprintf("%.2f",cx),family=font_family,size=1.9,colour="#1A1A1A")
  p }

# ---- manual horizontal legends (full width, so panels stay aligned) ----
make_hleg <- function(labs_v, cols, xs, xlim){
  df <- data.frame(lab=factor(labs_v, levels=labs_v), x=xs)
  ggplot(df) +
    geom_segment(aes(x=x-0.30, xend=x-0.08, y=1, yend=1, colour=lab), linewidth=0.5, lineend="round") +
    geom_point(aes(x=x-0.19, y=1, colour=lab), size=1.1) +
    geom_text(aes(x=x, y=1, label=lab), hjust=0, vjust=0.5, size=1.95, family=font_family, colour="black") +
    scale_colour_manual(values=cols, guide="none") +
    coord_cartesian(xlim=xlim, ylim=c(0.8,1.2), clip="off") +
    theme_void() + theme(plot.margin=margin(0,0,0,0,unit="mm")) }
leg_ind  <- make_hleg(c("AED","Et","ED","SMrz","SMs"), ind_colours,
                      xs=c(1.6,3.0,4.1,5.3,6.8), xlim=c(0.6,8.4))
# same xlim/scale as the top legend; "Atmospheric demand"->"AED" for consistency
# with the top legend and to keep the three items grouped & centred
leg_ctrl <- make_hleg(c("Soil moisture","Joint (non-separable)","AED"),
                      unname(ctrl_colours[c("Soil moisture","Joint (non-separable)","Atmospheric demand")]),
                      xs=c(1.7,3.6,6.1), xlim=c(0.6,8.4))
xlab_g <- grid::textGrob("Aridity index (AI)", gp=grid::gpar(fontsize=base_size, fontfamily=font_family, col="black"))

# ---- flat assembly (pure | and / operators -> columns align across rows) ----
top <- ptop(1) | ptop(2) | ptop(3) | ptop(4)
bot <- pbot(1) | pbot(2) | pbot(3) | pbot(4)
fig <- top / wrap_elements(leg_ind) / bot / wrap_elements(full=xlab_g) / wrap_elements(leg_ctrl) +
  plot_layout(heights=c(1, 0.10, 1, 0.05, 0.10)) &
  theme(plot.background=element_rect(fill="white",colour=NA))
ggsave(out_tif,fig,width=fig_width_mm,height=fig_height_mm,units="mm",dpi=fig_dpi,bg="white",device="tiff",compression="lzw")
#ggsave(out_pdf,fig,width=fig_width_mm,height=fig_height_mm,units="mm",dpi=fig_dpi,bg="white",device=cairo_pdf)
message("saved:\n  ", out_tif)
