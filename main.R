suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(CATALYST)
})

ctx = tercenCtx()

data <- ctx$as.matrix() %>% t()
colnames(data) <- ctx$rselect()[[1]]

files <- ctx$cselect() %>% 
  select(contains("filename"))

fset <- data %>% as_tibble() %>% bind_cols(files) %>%
  group_by({if("filename" %in% names(.)) filename else NULL}) %>% 
  group_map(~tim::matrix_to_flowFrame(as.matrix(.x))) %>%
  flowCore::flowSet()

sce <- prepData(fset, transform = FALSE)

res <- try(ctx$labels)

if (class(res) == "try-error"){ #specify mass channels stained for & debarcode
  channels <- ctx$op.value("channels", as.character, "")
  bc_ms <- eval(parse(text = paste0("c(", channels, ")")))
  if(is.null(bc_ms)) bc_ms <- as.numeric(gsub("[^0-9.-]", "", colnames(data)))
  
  sce <- assignPrelim(sce, bc_ms,assay = "counts", verbose = FALSE)
  sce <- estCutoffs(sce)
  sce <- applyCutoffs(sce,assay = "counts")
  
  # compute & extract spillover matrix
  sce <- computeSpillmat(
    sce,
    method = "default",
    trim = 0.5,
    th = 1e-05
  )
  
  sm <- metadata(sce)$spillover_matrix
} else {
  # get comp
  doc.id <- ctx$select(ctx$labels[[1]], nr = 1)[[1]]
  sm_df <- ctx$client$tableSchemaService$select(doc.id) %>%
    as_tibble()%>%as.data.frame()
  
  rownames(sm_df) <- sm_df[,1]
  sm_df[,1] <- NULL
  sm <-sm_df
  
  
  rowData(sce)$is_bc<-rowData(sce)$use_channel
  
  metadata(sce)$spillover_matrix<-as.matrix(sm)
}

sm<-as.matrix(sm)
sm[sm>1]<-1
p <- plotSpillmat(sce, sm = sm)

p_file <- suppressWarnings({tim::save_plot(p)})
df_plot <- tim::plot_file_to_df(p_file) %>%
  ctx$addNamespace() %>%
  as_relation()

sce <- compCytof(sce,sm = sm, method ="nnls" , assay = "counts", overwrite = FALSE, transform = FALSE, cofactor = 5)#"nnls" or "flow"
df <- assay(sce, "compcounts")

rids <- ctx$rselect()[1]

colnames(rids) <- "channel"

df_out <- df %>%
  as_tibble(rownames = "channel") %>%
  tidyr::pivot_longer(cols = !contains("channel"), names_to = ".ci") %>%
  mutate(.ci = as.integer(gsub("V", "", .ci)) - 1L) %>%
  left_join(rids %>% mutate(.ri = seq(1, nrow(.)) - 1L), by = "channel") %>%
  select(-channel) %>%
  ctx$addNamespace() %>%
  as_relation()


join_res = df_out %>%
  left_join_relation(ctx$crelation, ".ci", ctx$crelation$rids) %>%
  left_join_relation(df_plot, list(), list()) %>%
  as_join_operator(ctx$cnames, ctx$cnames)

 join_res %>%
  save_relation(ctx)
