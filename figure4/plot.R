# Figure 4 a-e; panel e raised symbols, no x/y grids, EGFR nominal-P footer.
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(grid); library(svglite); library(jsonlite)
})
OUT <- output_path(cfg, "figure4")
for (d in c("", "panels", "source_data", "qa"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
input_rds <- file.path(OUT, "rds/figure4_objects.rds")
before <- tools::md5sum(input_rds)
o <- readRDS(input_rds)

# Exact palette and axis styling from original Figure 4 build v6.
LINES <- c("MDA-MB-231", "BT-549", "MDA-MB-468")
PAL <- setNames(c("#4276B9", "#007C91", "#C57829"), LINES)
INK <- "#202E38"; MUTED <- "#667582"; RULE <- "#DAE1E5"
white <- function(z, a = .22) {
  v <- col2rgb(z) / 255
  rgb(1-a+a*v[1, ], 1-a+a*v[2, ], 1-a+a*v[3, ])
}
symbol <- function(q) ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", "ns")))
short <- function(z) unname(as.character(z)) # full cell-line names for display and exports
save_data <- function(d, n) write.csv(d, file.path(OUT, "source_data", paste0(n, ".csv")), row.names = FALSE)
theme_set(theme_classic(base_family = cfg$font_family, base_size = 8) +
  theme(text = element_text(face = "plain", colour = INK),
    axis.text = element_text(size = 7, colour = INK),
    axis.title = element_text(size = 7.5, face = "plain"),
    axis.line = element_line(colour = "#93A2AD", linewidth = .25),
    axis.ticks = element_line(colour = "#93A2AD", linewidth = .22),
    axis.ticks.x = element_blank(), axis.ticks.length = unit(.8, "mm"),
    panel.grid.major.y = element_line(colour = "#EAF0F3", linewidth = .18),
    panel.grid.minor = element_blank(), legend.position = "none",
    plot.title = element_text(size = 8, face = "italic", margin = margin(b = 3)),
    plot.margin = margin(1, 1, 1, 1, "mm")))

BASE <- c("OSMR", "IL6ST", "CDKN1A", "EGFR")
meta_a <- o$meta_acute %>% mutate(LineS = short(o$LINE[CellLine]), Model = "Acute", Cond = Condition)
meta_o <- o$meta_oe %>% mutate(LineS = short(o$LINE[CellLine]), Model = "OE", Cond = ifelse(Treatment == "OSM", "OSM", "V"))
long <- function(mat, meta, genes) {
  as.data.frame(t(mat[genes, meta$Sample, drop = FALSE])) %>%
    rownames_to_column("Sample") %>% pivot_longer(-Sample, names_to = "Gene", values_to = "Value") %>%
    left_join(meta, by = "Sample") %>%
    mutate(LineS = factor(LineS, levels = LINES))
}
x <- bind_rows(long(o$log_acute, meta_a, BASE),
               long(o$log_oe, meta_o, "CDKN1A"))
ci <- o$ci_all %>% filter(!grepl("Interaction", Contrast)) %>%
  mutate(LineS = factor(short(o$LINE[sub("_.*", "", Contrast)]), levels = LINES),
    Model = ifelse(grepl("Response", Contrast), "OE", "Acute"),
    Cond = ifelse(Model == "OE", "OSM", sub(".*_(\\d+h)_vs_NT", "\\1", Contrast)),
    Symbol = symbol(FDR)) %>% filter(!is.na(LineS))
stopifnot(n_distinct(ci$Contrast) == 8)

# Recalculate the historical bandwidths from TPM observations rather than
# depending on a cached CSV. Pool within-group residuals over acute AND
# sustained libraries for each displayed gene (original v4 rule).
bw <- bind_rows(long(o$log_acute, meta_a, BASE), long(o$log_oe, meta_o, BASE)) %>%
  group_by(Gene, Model, LineS, Cond) %>% mutate(resid = Value - mean(Value)) %>%
  ungroup() %>% group_by(Gene) %>%
  summarise(Bandwidth = max(.055, 1.06 * sd(resid) * 3^(-.2)), .groups = "drop")
save_data(bw, "violin_bandwidths")
save_data(x, "expression_observations")
point_log <- list(); density_log <- list(); summary_log <- list()

groups <- data.frame(Model = c(rep("Acute", 3), rep("OE", 2)),
  LineS = c(LINES, LINES[c(1, 3)]), Offset = c(0, 3.5, 7, 11, 13.7), N = c(3, 3, 3, 2, 2))
groups$Mid <- groups$Offset + (groups$N + 1)/2
groups$Colour <- PAL[groups$LineS]
save_data(groups, "condition_positions")

base_stats <- o$bl_stat %>% filter(Gene %in% BASE) %>%
  mutate(x1 = pmin(match(short(A), LINES), match(short(B), LINES)),
         x2 = pmax(match(short(A), LINES), match(short(B), LINES)),
         Symbol = symbol(p_holm), Row = ifelse(x1 == 1 & x2 == 2, 1, ifelse(x1 == 2 & x2 == 3, 2, 3)))
save_data(base_stats, "a_baseline_Welch_Holm_from_B3")
save_data(o$bl_stat, "baseline_complete_18_test_family")
save_data(filter(x, Gene %in% BASE, Model == "Acute", Cond == "NT"), "a_baseline_observations")

make_violin <- function(gene, mode = c("baseline", "response"), y_max, show_y = TRUE) {
  mode <- match.arg(mode)
  d <- filter(x, Gene == gene)
  if (mode == "baseline") d <- filter(d, Model == "Acute", Cond == "NT")
  d <- d %>% left_join(groups, by = c("Model", "LineS")) %>%
    group_by(Model, LineS, Cond) %>% arrange(Sample, .by_group = TRUE) %>%
    mutate(cx = if (mode == "baseline") match(LineS, LINES) else
      Offset + match(Cond, ifelse(Model[1] == "OE", list(c("V", "OSM")), list(c("NT", "6h", "24h")))[[1]]),
      px = cx + seq(-.09, .09, length.out = n()),
      pch = if (mode == "baseline") match(LineS, LINES)+20 else ifelse(Cond %in% c("NT", "V"), 21, ifelse(Cond == "6h", 22, 24)),
      Fill = ifelse(Cond %in% c("NT", "V"), "white", Colour)) %>% ungroup()
  sm <- d %>% group_by(Model, LineS, Cond, cx) %>%
    summarise(Mean = mean(Value), SD = sd(Value), n = n(), .groups = "drop")
  b <- bw$Bandwidth[bw$Gene == gene]
  den <- bind_rows(lapply(split(d, interaction(d$Model, d$LineS, d$Cond, drop = TRUE)), function(z) {
    de <- density(z$Value, bw = b, from = max(0, min(z$Value) - 3*b), to = max(z$Value) + 3*b, n = 160)
    wi <- de$y/max(de$y)*.28
    data.frame(X = c(z$cx[1]-wi, rev(z$cx[1]+wi)), Y = c(de$x, rev(de$x)),
      Group = paste(z$Model[1], z$LineS[1], z$Cond[1]), Colour = z$Colour[1], Fill = white(z$Colour[1]))
  }))
  key <- paste(mode, gene, sep = "_")
  point_log[[key]] <<- mutate(d, PlotID = key)
  density_log[[key]] <<- mutate(den, PlotID = key)
  summary_log[[key]] <<- mutate(sm, Gene = gene, PlotID = key)
  tk <- d %>% distinct(cx, Cond, LineS) %>% arrange(cx)
  labels <- if (mode == "baseline") as.character(tk$LineS) else c(NT="NT", `6h`="6 h", `24h`="24 h", V="FG", OSM="FG/OSM")[tk$Cond]
  xmax <- if(mode == "baseline") 3.5 else 16.2
  y_breaks <- pretty(c(0,y_max),n=4)
  p <- ggplot() +
    geom_polygon(data=den, aes(X,Y,group=Group,colour=Colour,fill=Fill), linewidth=.28) +
    geom_errorbar(data=sm, aes(cx,ymin=Mean-SD,ymax=Mean+SD), width=.13, linewidth=.28, colour=INK) +
    geom_segment(data=sm, aes(x=cx-.17,xend=cx+.17,y=Mean,yend=Mean), linewidth=.4, colour=INK) +
    geom_point(data=d, aes(px,Value,colour=Colour,fill=Fill,shape=pch), size=1.5, stroke=.34) +
    scale_colour_identity() + scale_fill_identity() + scale_shape_identity() +
    scale_x_continuous(breaks=tk$cx,labels=labels,limits=c(.5,xmax),expand=c(0,0)) +
    scale_y_continuous(breaks=y_breaks,expand=c(0,0)) +
    coord_cartesian(ylim=c(-.15,y_max),clip="off") +
    labs(x=NULL,y=if(show_y)expression(log[2](TPM+1))else NULL)
  if(mode == "baseline") {
    ann <- filter(base_stats,Gene==gene) %>% mutate(yy=y_max*c(.72,.82,.92)[Row])
    p <- p + geom_segment(data=ann,aes(x=x1,xend=x2,y=yy,yend=yy),linewidth=.25,colour=INK) +
      geom_segment(data=ann,aes(x=x1,xend=x1,y=yy-.12,yend=yy),linewidth=.25,colour=INK) +
      geom_segment(data=ann,aes(x=x2,xend=x2,y=yy-.12,yend=yy),linewidth=.25,colour=INK) +
      geom_text(data=ann,aes(x=(x1+x2)/2,y=yy+.015*y_max,label=Symbol),size=2.4,vjust=0,colour=INK) +
      labs(title=gene) + theme(axis.text.x=element_text(angle=55,hjust=1,vjust=1,size=6.5))
  } else {
    at <- ci %>% filter(Gene==gene) %>% inner_join(groups,by=c("Model","LineS"))
    at <- at %>% mutate(cx=Offset+ifelse(Model=="OE",2,match(Cond,c("NT","6h","24h"))))
    gh <- groups
    p <- p + geom_text(data=at,aes(cx,y_max*.91,label=Symbol),size=2.5,colour=INK) +
      geom_segment(data=gh,aes(x=Offset+.65,xend=Offset+N+.35,y=y_max*1.055,yend=y_max*1.055,colour=Colour),linewidth=.25) +
      geom_text(data=gh,aes(Mid,y_max*1.11,label=LineS,colour=Colour),size=2.3) +
      theme(plot.margin=margin(7.5,1,1,1,"mm"))
    if(mode=="response") p <- p +
      annotate("text",x=5.5,y=y_max*1.24,label="Acute OSM",size=2.65,colour=INK) +
      annotate("text",x=13.6,y=y_max*1.24,label="OSM expression",size=2.65,colour=INK) +
      annotate("segment",x=10.85,xend=10.85,y=-.15,yend=y_max*1.16,colour=RULE,linewidth=.25) +
      theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=6.5),plot.margin=margin(13,1,1,1,"mm"))
  }
  p
}

a_plots <- setNames(lapply(BASE, function(g) make_violin(g, "baseline", switch(g, OSMR=11.5, IL6ST=11, CDKN1A=10, EGFR=15), show_y=(g=="OSMR"))), BASE)
d_plot <- make_violin("CDKN1A", "response", 10)

save_data(filter(ci,Gene=="CDKN1A"),"d_CDKN1A_statistics")


# Differential-expression breadth. Preserve the supplied threshold exactly.
deg <- ci %>% filter(FDR < .05,abs(logFC) > .585) %>%
  mutate(Time=ifelse(Model=="OE","OE",Cond),Direction=ifelse(logFC>0,"Up","Down")) %>%
  count(LineS,Time,Direction,name="n") %>%
  mutate(Time=factor(Time,levels=c("6h","24h","OE")),Signed_n=ifelse(Direction=="Down",-n,n))
stopifnot(nrow(deg)==16)
save_data(deg,"b_DEG_counts")
count_cols <- c(Up="#D4A58D",Down="#98B4CF")
b_plot <- ggplot(deg,aes(Time,Signed_n,fill=Direction)) +
  geom_col(width=.62,colour=INK,linewidth=.2) +
  geom_hline(yintercept=0,colour="#93A2AD",linewidth=.25) +
  geom_text(aes(label=n,vjust=ifelse(Direction=="Down",1.25,-.45)),size=2.1,colour=INK) +
  facet_wrap(~LineS,nrow=1,drop=FALSE) +
  scale_fill_manual(values=count_cols) +
  scale_y_continuous(breaks=c(-1000,0,1000),labels=abs,limits=c(-1550,1650),expand=c(0,0)) +
  scale_x_discrete(drop=FALSE,labels=c(`6h`="6 h",`24h`="24 h",OE="FG/OSM")) +
  geom_text(data=data.frame(LineS=factor("BT-549",levels=LINES),Time=factor("OE",levels=c("6h","24h","OE"))),
    aes(x=Time,y=0,label="NA"),inherit.aes=FALSE,size=2.3,colour=MUTED,vjust=-.5) +
  labs(x=NULL,y="DEGs") +
  theme(strip.background=element_blank(),strip.text=element_text(size=6.5,face="plain",margin=margin(b=3)),
        panel.spacing=unit(2,"mm"),axis.text.x=element_text(size=6.3,angle=45,hjust=1,vjust=1),plot.margin=margin(1,1,1,0,"mm"))

# The same eight rows and chronological order are used for c and e.
rows <- data.frame(Contrast=c("MDA231_6h_vs_NT","MDA231_24h_vs_NT","BT549_6h_vs_NT","BT549_24h_vs_NT",
  "MDA468_6h_vs_NT","MDA468_24h_vs_NT","MB231_Response","MB468_Response"),
  RowY=c(8.5,7.5,6.3,5.3,4.1,3.1,1.6,.6),
  Label=c("MDA-MB-231  6 h","MDA-MB-231  24 h","BT-549  6 h","BT-549  24 h","MDA-MB-468  6 h","MDA-MB-468  24 h","MDA-MB-231  FG/OSM","MDA-MB-468  FG/OSM"))
gs <- o$gsea %>% filter(Hallmark=="IL6_JAK_STAT3_SIGNALING") %>%
  inner_join(rows,by="Contrast") %>%
  mutate(LineS=factor(short(Line),levels=LINES),Symbol=symbol(padj))
stopifnot(nrow(gs)==8,all(gs$NES>0))
save_data(gs,"c_IL6_JAK_STAT3_GSEA")
c_plot <- ggplot(gs,aes(NES,RowY,colour=LineS)) +
  geom_segment(aes(x=0,xend=NES,yend=RowY),linewidth=.45) +
  geom_point(aes(shape=LineS),size=1.8,stroke=.3) +
  geom_text(aes(x=NES+.12,label=Symbol),size=2.3,hjust=0,colour=INK) +
  geom_hline(yintercept=2.35,colour=RULE,linewidth=.25,linetype="dashed") +
  scale_colour_manual(values=PAL) + scale_shape_manual(values=c(16,15,17)) +
  scale_x_continuous(limits=c(0,3.2),breaks=0:3,expand=c(0,0)) +
  scale_y_continuous(limits=c(0,9.15),breaks=rows$RowY,labels=rows$Label,expand=c(0,0)) +
  labs(x="NES",y=NULL) +
  theme(axis.text.y=element_text(size=6.8),axis.ticks.y=element_blank(),
    panel.grid.major.y=element_blank(),plot.margin=margin(1,2,1,1,"mm"))

fo <- ci %>% filter(Gene %in% c("CDKN1A","EGFR")) %>% inner_join(rows,by="Contrast")
stopifnot(nrow(fo)==16,all(fo$lo<=fo$logFC),all(fo$hi>=fo$logFC))
save_data(fo,"e_effects_95CI_FDR")
fo <- fo %>% mutate(Pos = match(Contrast, rows$Contrast),
  XLabel = paste(LineS, ifelse(Model == "OE", "FG/OSM", ifelse(Cond == "6h", "6 h", "24 h"))))
save_data(fo, "e_effects_95CI_FDR_P")
save_data(filter(fo, Gene == "EGFR") %>% arrange(Pos), "e_EGFR_nominal_P_footer")
effect_plot <- function(g) {
  d <- filter(fo, Gene == g) %>% arrange(Pos)
  ggplot(d, aes(Pos, logFC, colour = LineS)) +
    geom_hline(yintercept = 0, colour = "#93A2AD", linewidth = .30) +
    geom_vline(xintercept = 6.5, colour = RULE, linewidth = .25, linetype = "dashed") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .21, linewidth = .35) +
    geom_point(aes(shape = LineS), size = 1.9, stroke = .3) +
    geom_text(aes(y = hi + .26, label = Symbol), size = 2.35, vjust = 0, colour = INK) +
    scale_colour_manual(values = PAL) + scale_shape_manual(values = c(16, 15, 17)) +
    scale_x_continuous(limits = c(.55, 8.5), breaks = 1:8, labels = d$XLabel, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.12, 3.5), breaks = -1:3, expand = c(0, 0)) +
    labs(title = g, x = NULL, y = expression(Delta~log[2](TPM+1))) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 6.5),
      axis.text.y = element_text(size = 7),
      panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
      plot.title = element_text(size = 8.5, face = "italic", hjust = .5, margin = margin(b = 4)),
      plot.margin = margin(1, 3, 1, 1, "mm"))
}
e_plots <- list(effect_plot("CDKN1A"), effect_plot("EGFR"))
save_data(fo %>% transmute(Gene, Contrast, Symbol, CI_upper = hi, Label_bottom = hi + .26), "e_annotation_positions")

# Fixed physical layout makes panel exports and the assembled figure identical.
txt <- function(s,x,y,size=8,col=INK,just="left",font="plain")
  grid.text(s,x=unit(x,"mm"),y=unit(1,"npc")-unit(y,"mm"),just=c(just,"top"),
    gp=gpar(fontfamily=cfg$font_family,fontsize=size,fontface=font,col=col))
place <- function(p,x,y,w,h) {
  pushViewport(viewport(x=unit(x,"mm"),y=unit(1,"npc")-unit(y,"mm"),width=unit(w,"mm"),height=unit(h,"mm"),just=c("left","top")))
  grid.draw(ggplotGrob(p));popViewport()
}
head <- function(letter,title) {txt(paste0(letter,")"),0,0,10.5,font="bold");txt(title,5.8,.55,8.1)}
draw_a <- function() {
  head("a","Baseline expression")
  for(i in seq_along(BASE))place(a_plots[[i]],c(0,28,54,80)[i],7,c(28,26,26,26)[i],49)
}
draw_b <- function() {
  head("b","OSM-responsive genes")
  for(i in 1:2) {
    xx<-c(21,40)[i];lab<-c("Up","Down")[i]
    grid.rect(x=unit(xx,"mm"),y=unit(1,"npc")-unit(7.5,"mm"),width=unit(2.4,"mm"),height=unit(2,"mm"),gp=gpar(fill=count_cols[lab],col=INK,lwd=.4))
    txt(lab,xx+2,6.3,6.6)
  }
  place(b_plot,0,11,68,44)
}
draw_c <- function() {
  head("c","IL6 / JAK / STAT3")
  place(c_plot,0,6.5,63,47)
}
draw_d <- function() {
  head("d","");txt("CDKN1A",5.8,.55,8.1,font="italic")
  place(d_plot,0,6.5,111,47)
}
draw_e <- function() {
  head("e", "Transcript effect estimates")
  place(e_plots[[1]], 0, 7, 88, 56)
  place(e_plots[[2]], 90, 7, 88, 56)
  txt("Pointwise 95% CI (unadjusted); symbols: genome-wide FDR", 2, 65.5, 6.7, col = MUTED)
  txt("EGFR unadjusted P (two-sided moderated t tests):", 2, 70.0, 6.7)
  ep <- filter(fo, Gene == "EGFR")
  for(i in seq_along(LINES)) {
    vals <- ep %>% filter(LineS == LINES[i]) %>% arrange(Pos)
    labs <- ifelse(vals$Model == "OE", "FG/OSM", ifelse(vals$Cond == "6h", "6 h", "24 h"))
    v <- ifelse(vals$P < .1, sprintf("%.4f", vals$P), sprintf("%.3f", vals$P))
    label <- paste0(LINES[i], ": ", paste(paste0(labs, " = ", v), collapse = "; "))
    txt(label, 2, 74.0 + (i - 1)*3.4, 6.5, col = MUTED)
  }
}
layout <- data.frame(Panel=letters[1:5],X=c(2,110,2,67,2),Y=c(2,2,59,59,119),W=c(106,68,63,111,178),H=c(56,56,54,54,86))
drawers <- list(draw_a,draw_b,draw_c,draw_d,draw_e)
save_data(layout, "layout_mm")
W <- 182; H <- 207
draw_main <- function() {
  grid.newpage()
  for(i in 1:5) {
    z <- layout[i,]
    pushViewport(viewport(x=unit(z$X,"mm"),y=unit(1,"npc")-unit(z$Y,"mm"),width=unit(z$W,"mm"),height=unit(z$H,"mm"),just=c("left","top")))
    drawers[[i]]();popViewport()
  }
}
export <- function(draw,base,w,h,dpi=600) {
  cairo_pdf(paste0(base,".pdf"),width=w/25.4,height=h/25.4,family=cfg$font_family,bg="white")
  draw();dev.off()
  svglite(paste0(base,".svg"),width=w/25.4,height=h/25.4,bg="white",system_fonts=list(sans=cfg$font_family))
  draw();dev.off()
  png(paste0(base,".png"),width=w,height=h,units="mm",res=dpi,type="cairo",bg="white")
  draw();dev.off()
}
export(draw_main,file.path(OUT,"Figure4_a-e"),W,H)
for(i in 1:5)export(function(){grid.newpage();drawers[[i]]()},file.path(OUT,"panels",paste0("Figure4",letters[i])),layout$W[i],layout$H[i],400)
save_data(bind_rows(point_log),"plotted_points")
save_data(bind_rows(density_log),"violin_vertices")
save_data(bind_rows(summary_log),"plotted_mean_SD")
save_data(data.frame(LineS=LINES,Colour=PAL),"original_figure4_palette")
after<-tools::md5sum(input_rds)
stopifnot(identical(before,after))
write.csv(data.frame(Path=names(before),MD5=unname(before)),file.path(OUT,"qa/input_md5.csv"),row.names=FALSE)
write_json(list(status="PASS",panels=letters[1:5],excluded=c("f","g"),dimensions_mm=c(W,H),
  input_files_unchanged=TRUE,contrasts=8,CI_rows=nrow(fo),GSEA_rows=nrow(gs),DEG_count_rows=nrow(deg),
  sample_sizes=x %>% distinct(Sample,Model,LineS,Cond) %>% count(Model,LineS,Cond),
  EGFR_FDR_all_ge_005=all(filter(fo,Gene=="EGFR")$FDR>=.05),
  EGFR_CI_excluding_zero=filter(fo,Gene=="EGFR",lo>0|hi<0) %>% select(Contrast,lo,hi,FDR)),
  file.path(OUT,"qa/validation.json"),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(OUT,"qa/sessionInfo.txt"))
cat("Figure 4 a-e -> ",OUT,"\n",sep="")
