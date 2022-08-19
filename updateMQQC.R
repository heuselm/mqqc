#' ---
#' title: "QC result plotting from MaxQuant results"
#' author: mhe
#' output:
#'   html_document:
#'     keep_tex: false
#' ---

#'## Environment setup
#+ r load packages
library('data.table')
library('ggplot2')
library('plotly')
library('ggrepel')
library('htmlwidgets')

# evosep color scheme
escolors = c("#87B09A", "#FF671B", "#D1E0D7", "#222720", "#8A8A8D")
n_latest = 5

# get working directory name as name of the MS instrument
MSname = basename(getwd())

# copy/backup current plots to history
for (file in list.files(pattern = "^QC_")){
  print(file)
  print(paste0("To: ", "QC_plots_history/", Sys.Date(), file))
  file.copy(from = file, paste0("QC_plots_history/", Sys.Date(), "_", file))
}
file.copy(from = "runAnnotation.csv", paste0("QC_plots_history/", Sys.Date(), "_runAnnotation.csv"))


#'## Import data
#+ r Import data
# load all /mq_txt_results/txt_
dirs = list.dirs(path = ".//mq_txt_results")
ndirs = length(list.dirs(path = ".//mq_txt_results"))
result_dirs = dirs[2:ndirs]

ev_combined = data.table()
summary_combined = data.table()

for (resultset in result_dirs){
  ## collect evidence files
  ev = fread(paste0(resultset, "/evidence.txt"))
  # clean column names
  names(ev) = gsub("\\]|\\[|\\)|\\(|\\/| ", "_", names(ev))
  names(ev) = gsub("_$", "", names(ev))
  names(ev) = gsub("__", "_", names(ev))
  names(ev) = gsub("%", "percent", names(ev))
  ev[, analysis_directory:=basename(resultset)]
  ev_combined = rbind(ev_combined, ev, fill = TRUE)
  
  ## collect summary files
  summ = fread(paste0(resultset, "/summary.txt"))
  names(summ) = gsub("\\]|\\[|\\)|\\(|\\/| ", "_", names(summ))
  names(summ) = gsub("_$", "", names(summ))
  names(summ) = gsub("__", "_", names(summ))
  names(summ) = gsub("%", "percent", names(summ))
  
  summ[, analysis_directory:=basename(resultset)]
  summary_combined = rbind(summary_combined, summ, fill = TRUE)
}

# write a template runAnnotation.txt table
ev_combined[, Acq_day:=as.numeric(unlist(strsplit(analysis_directory, split = "_"))[2]), Raw_file]
annotationTemplate = unique(ev_combined[, .(Raw_file, Acq_day)])
annotationTemplate[, c("operator", "LC_id",	"LC_column",	"LC_meth", "Sample_id", "Sample_load", "QC_boxid_loaddate") := list("FILL_IN")]
fwrite(annotationTemplate, file = "runAnnotation_template.csv")
ev_combined$Acq_day = NULL

# read & merge current runAnnotation
ann = fread("runAnnotation.csv")
nruns = length(unique(ev_combined$Raw_file))
ev_combined = merge(ev_combined, ann, by = "Raw_file")

# define n latest n QC analyses:
latest_n = sort(unique(ev_combined$Acq_day), decreasing = T)[1:n_latest]

# View(unique(ev_combined[, .(Raw_file, Acq_day, analysis_directory)]))

#'## Count & plot IDs
#+ r Count & plot IDs
ids.prot = ev_combined[, length(unique(Proteins)), by = names(ann)]
ids.pep = ev_combined[, length(unique(Modified_sequence)), by = names(ann)]
ids = rbind(ids.prot[, level:="Protein.groups"],
                  ids.pep[, level:="Peptides"])

names(ids)

p1 = ggplot(ids,
            aes(x = as.character(Acq_day), y = V1, col = Sample_id, group = Sample_id)) +
  geom_point() +
  geom_line() +
  ggtitle(paste(MSname, "QC Identifications MaxQuant")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  #geom_text(aes(x = as.character(Acq_day), y = V1, label=V1), vjust = -2, cex = 4)  +
  ylab("N") +
  facet_wrap(level~LC_meth, scales = "free") +
  theme(axis.text.x=element_text()) + 
  expand_limits(y=0)
p1 
ggsave("QC_IDs_lineplots.pdf", width = (nruns/8)+5, limitsize = F)
dev.off()
htmlwidgets::saveWidget(ggplotly(p1), file = "QC_IDs_lineplots.html")
# ggplotly(p1)

#'## Mass errors
#+ r Mass errors
names(ev_combined)
p2 = ggplot(ev_combined, aes(text = Raw_file)) +
  geom_boxplot(aes(x = as.character(Acq_day), y = Mass_error_ppm, color = Sample_id, group = Raw_file)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(LC_meth~Acq_day,  scales = "free_x") +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("Mass errors [ppm]")
p2
ggsave("QC_Mass_error_ppm.pdf", height = length(unique(ev_combined$LC_meth))*3, width = (ndirs/4)+3)
htmlwidgets::saveWidget(ggplotly(p2), file = "QC_Mass_error_ppm.html")
# ggplotly(p2)

#'## ID over RT density plots / 'IDTic'
#+ r ID over RT density plots
names(ev_combined)
p3 = ggplot(ev_combined[Acq_day %in% latest_n]) +
  geom_density(aes(x = Retention_time, y = ..count.., color = as.character(Acq_day), group = Raw_file), adjust = 0.2) +
  theme_minimal() + 
  facet_grid(~LC_meth, scales = "free") +
  ggtitle(paste(n_latest, "latest QC sets - Identification density along RT"))
p3
ggsave("QC_ID_density_along_RT.pdf", height = length(unique(ev_combined$LC_meth))*3, width = (ndirs/4)+3)
htmlwidgets::saveWidget(ggplotly(p3), file = "QC_ID_density_along_RT.html")
# ggplotly(p3)

#'## Peak shape / retention length per RT segment
#+ Peak shape / retention length per RT segment
names(ev_combined)
p4 = ggplot(ev_combined[Acq_day %in% latest_n]) +
  geom_boxplot(aes(x = cut(Retention_time, 5), y = Retention_length, color = Raw_file),
               position = position_dodge2(preserve = "single"),
               outlier.shape = NA) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(LC_meth~Acq_day,  scales = "free") +
  theme(legend.position = "none") +
  ggtitle("Peak width per RT section") +
  xlab("Retention time range")
p4
ggsave("QC_Retention_length.pdf", height = length(unique(ev_combined$LC_meth))*3, width = (nruns/4)+3)
htmlwidgets::saveWidget(ggplotly(p3), file = "QC_Retention_length.html")
# ggplotly(p4)

#'## RT of specific HeLa peptides DSYVGDEAQSK_599.765++, GYSFTTTAER_566.76++, TVTAMDVVYALK_655.855++
#+ r RT of specific HeLa peptides
targets_hela = c("DSYVGDEAQSK", "GYSFTTTAER", "TVTAMDVVYALK")
p5 = ggplot(ev_combined[Sequence %in% targets_hela],
            aes(x = as.character(Acq_day),
                y = Retention_time,
                col = paste(Modified_sequence,
                            Charge, m_z),
                group = paste(Modified_sequence, Charge),
                text = Raw_file)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(LC_id~LC_meth) +
  ggtitle(paste(MSname, "3 peptides' retention times in QC runs"))
p5
ggsave("QC_HeLa_3peptides_RTs.pdf", height = 8, width = 10)
htmlwidgets::saveWidget(ggplotly(p5), file = "QC_HeLa_3peptides_RTs.html")
# ggplotly(p5)


#'## MaxQuant Summary
#+ r MaxQuant Summary
names(summary_combined)
# merge annotation
summary_combined_ann = merge(summary_combined, ann, by = "Raw_file")
# write combined csv table
fwrite(summary_combined_ann, file = "QC_Summary_table.csv")
# select columns for graphs
names(summary_combined_ann)
selected_columns = c(names(summary_combined)[grep("^MS|Peak|Peptide|deviation", names(summary_combined))])
# plot graphs
pdf("QC_Summary_graphs.pdf", height = 5, width = 3+length(unique(summary_combined_ann$Acq_day))/4)
for (col in selected_columns) {
  print(col)
  p = ggplot(summary_combined_ann, aes_string(x = "as.character(Acq_day)", y = col, col = "LC_id")) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_wrap(~LC_meth) +
    coord_cartesian(ylim = c(0, 1.1*max(summary_combined_ann[[col]]))) +
    ggtitle(paste("MaxQuant Summary: ", col, "@", MSname))
  plot(p)
}
dev.off()
plot.new()

