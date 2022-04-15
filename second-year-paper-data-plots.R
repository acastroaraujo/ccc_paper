

# Setup -------------------------------------------------------------------

library(tidyverse)
library(igraph)
library(ggraph)
library(backbone)
library(stm)
library(loo)
library(patchwork)

theme_set(
  theme_light(base_family = "Amiri") + theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(colour = "white"), 
    strip.background = element_rect(fill = "#4C4C4C")
  )
)


# Data --------------------------------------------------------------------

get_data <- function(x) {
  stopifnot(x %in% c("citations", "descriptors", "metadata"))
  path <- "~/Documents/Repositories/ccc_datos/data/"
  readr::read_rds(stringr::str_glue("{path}{x}.rds"))
}

el <- get_data("citations")

metadata <- get_data("metadata") |> 
  select(id, type, year, date)

# Clean up

metadata <- metadata |> 
  filter(!id %in% c("T-680-07", "T-277-09"))

el <- el |> 
  filter(from %in% metadata$id, to %in% metadata$id) |> 
  filter(from_date > to_date) |>       ## remove time travel (this can happen in a few cases)
  distinct(from, to, .keep_all = TRUE) ## remove multiple citations per case

dag <- igraph::graph_from_data_frame(
  d = el, 
  directed = TRUE,
  vertices = metadata
)

metadata <- metadata |> 
  left_join(
    igraph::degree(dag, mode = "in") |> 
      tibble::enframe("id", "indegree")
  ) |> 
  left_join(
    igraph::degree(dag, mode = "out") |> 
      tibble::enframe("id", "outdegree")
  ) 

cd_index <- read_rds(here::here("analysis", "out", "cd_index.rds"))
atypicality <- read_rds(here::here("analysis", "out", "atypicality.rds"))

adumbration <- read_rds(here::here("analysis", "out", "adumbration.rds"))
adumbration_counts <- read_rds(here::here("analysis", "out", "adumbration_counts.rds"))
df_threshold <- adumbration_counts |> 
  mutate(year = ccc:::extract_year(id)) |> 
  group_by(year) |> 
  summarize(threshold = unique(threshold))

transposition <- read_rds(here::here("analysis", "out", "transposition.rds"))

metadata <- metadata |> 
  left_join(atypicality) |> 
  left_join(adumbration_counts |> select(!threshold))

metadata <- metadata |> 
  mutate(across(c(adum_count, adum_max), \(x) ifelse(is.na(x), 0, x)))

metadata <- metadata |> 
  left_join(transposition)

metadata <- metadata |> 
  mutate(days_of_exposure = as.integer(max(date) - date))

mod_coeffplot <- read_rds(here::here("analysis", "out", "mod_coeffplot.rds"))
loo_out <- read_rds(here::here("analysis", "out", "loo_out2.rds"))
mod8_ppc_plots <- read_rds(here::here("analysis", "out", "mod8_ppc_plots.rds"))

corr_out <- read_rds(here::here("analysis", "out", "corr_tm.rds"))

## recombination

bb <- read_rds(here::here("analysis", "out", "bb_sdsm_principles.rds"))
ALPHA <- 0.0001
bb_net <- backbone.extract(bb, alpha = ALPHA, class = "igraph", signed = TRUE)

# bb <- read_rds(here::here("analysis", "out", "bb_rp_sdsm.rds"))
# bb_net <- igraph::graph_from_adjacency_matrix(bb, mode = "undirected", diag = FALSE, weighted = TRUE)




bb_net_positive <- delete_edges(bb_net, which(E(bb_net)$weight == -1)) 
#bb_net_positive <- delete_vertices(bb_net_positive, which(degree(bb_net_positive) == 0))
igraph::V(bb_net_positive)$degree <- igraph::degree(bb_net_positive)
igraph::V(bb_net_positive)$b <- igraph::betweenness(bb_net_positive, directed = FALSE)

igraph::V(bb_net_positive)$cluster <- igraph::cluster_louvain(bb_net_positive) |> 
  igraph::membership() |> 
  factor()

principles_target <- c(
  "principio de laicidad", "principio de pluralismo religioso",
  #"principio de certeza tributaria", "principio de irretroactividad tributaria",
  #"principio unitario del estado","principio de autonomia territorial",
  #"principio de unidad de materia", "principio de publicidad",
  "principio de inmediatez", "principio de subsidiariedad",
  "principio pro infans", "principio de corresponsabilidad"
)

V(bb_net)$cluster <- V(bb_net_positive)$cluster
V(bb_net)$degree <- V(bb_net_positive)$degree
V(bb_net)$b <- V(bb_net_positive)$b

pp_mat <- read_rds(here::here("analysis", "out", "principles_mat.rds"))
#pmat <- read_rds(here::here("analysis", "out", "principles_mat.rds"))
#rmat <- read_rds(here::here("analysis", "out", "rights_mat.rds"))
#pp_mat <- cbind(pmat, rmat)
case_principles_keep <- names(pp_mat["C-370-06", ][which(pp_mat["C-370-06", ] > 0)])

## transposition 

stm_corr <- cor(corr_out$theta) ## topic correlations
colnames(stm_corr) <- rownames(stm_corr) <- 1:nrow(stm_corr)

stm_net <- igraph::graph_from_adjacency_matrix(
  adjmatrix = stm_corr, 
  weighted = TRUE, 
  mode = "undirected", 
  diag = FALSE
)

stm_el <- igraph::as_data_frame(stm_net) |> 
  rename(corr = weight) ## corelation edge list

corr_net <- stminsights::get_network(corr_out, cutoff = 0.15, labels = 1:ncol(corr_out$theta)) 
V(corr_net)$d <- colSums(corr_out$theta)
#corr_net <- delete.vertices(corr_net, which(degree(corr_net) == 0))
V(corr_net)$strength <- strength(corr_net)
V(corr_net)$degree <- degree(corr_net)
V(corr_net)$betweenness <- betweenness(corr_net)

igraph::V(corr_net)$cluster <- igraph::cluster_louvain(corr_net) |> 
  igraph::membership() |> 
  as.character()

ccc_sparse <- read_rds(here::here("analysis", "out", "sparse_dtm.rds"))

row_remove <- which(!rownames(ccc_sparse) %in% metadata$id)
ccc_sparse <- ccc_sparse[-row_remove, ]
rownames(corr_out$theta) <- rownames(ccc_sparse)

index <- which(rownames(corr_out$theta) == "T-859-03")

x <- corr_out$theta[index, , drop = FALSE]
colnames(x) <- seq_along(x)
ego <- t(x) %*% x

ego_el <- igraph::graph_from_adjacency_matrix(
  adjmatrix = ego, 
  mode = "undirected", weighted = TRUE, diag = FALSE
) |> 
  igraph::as_data_frame() |> 
  mutate(across(c(to, from), as.character))

case_transposition <- full_join(ego_el, stm_el) |> 
  mutate(score = weight*corr) |> 
  filter(
    score < quantile(score, probs = 0.001) | score > quantile(score, probs = 0.999)
  ) 




# Plots -------------------------------------------------------------------

avg_per_year_df <- metadata |> 
  count(year, type) |> 
  group_by(type) |> 
  summarize(avg = mean(n)) |> 
  mutate(offset = case_when(
    type == "C" ~ 20,
    type == "SU" ~ 2,
    type == "T" ~ 60
  )) 

avg_per_year_df$x <- c(2019, 2008, 2019)

## plots
figure1 <- metadata |> 
  count(year, type) |>
  group_by(type) |> 
  mutate(avg = mean(n)) |> 
  ggplot(aes(year, n)) +
  geom_col(width = 3/4, position = "dodge") +
  geom_label(
    data = avg_per_year_df,
    mapping = aes(x = x, y = avg + offset, label = scales::comma(avg, 0.1)),
    size = 3
  ) +
  geom_hline(
    data = avg_per_year_df,
    mapping = aes(yintercept = avg), 
    linetype = "dashed"
  ) + 
  facet_wrap(~type, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = NULL, title = "Number of decisions per year",
       subtitle = "1992-2021",
       caption = "dashed lines indicate average decisions per year")

figure2 <- metadata |> 
  pivot_longer(c(indegree, outdegree), names_to = "dtype", values_to = "degree") |> 
  mutate(dtype = case_when(
    dtype == "indegree" ~ "Average Inward Citations (in-degree)",
    dtype == "outdegree" ~ "Average Outward Citations (out-degree)"
  )) |> 
  ggplot(aes(year, degree)) + 
  stat_summary(
    fun.data = mean_cl_boot, fatten = 1, 
    fun.args = list(conf.int = 0.9)
  ) + 
  facet_grid(type~dtype, scales = "free_y") + 
  scale_x_continuous(labels = seq(1992, 2022, 4), breaks = seq(1992, 2022, 4)) +
  labs(y = NULL, x = NULL) +
  theme(strip.text.y = element_text(angle = 0))

figure3 <- tidyr::crossing(
  x = seq(0, 5, 1),
  year = seq(1992, 2021)
) |> 
  mutate(prop = map2_dbl(x, year, function(x, y) {
    mean(metadata$indegree[which(metadata$year == y)] <= x)
  })) |> 
  filter(x %in% c(0)) |> 
  ggplot(aes(year, prop)) + 
  geom_col() + 
  geom_text(aes(label = scales::percent(prop, 1)), size = 2, nudge_y = 0.025) +
  scale_x_continuous(labels = seq(1992, 2021, 2), breaks = seq(1992, 2021, 2)) +
  theme(axis.text.y = element_blank()) +
  labs(
    x = NULL, y = NULL,
    title = "Proportion of cases that have never been cited",
    subtitle = paste0("As of ", max(metadata$date))
  )

figure4 <- adumbration |> 
  mutate(year = ccc:::extract_year(target)) |>
  left_join(df_threshold) |> 
  group_by(target) |> 
  filter(target %in% c("C-776-03", "C-370-06", "T-444-13")) |> 
  mutate(r_label = ifelse(days > threshold, reference, NA_character_)) |> 
  ggplot(aes(date, days)) + 
  geom_point(aes(color = target), show.legend = FALSE) + 
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  ggrepel::geom_text_repel(aes(label = r_label), size = 3, direction = "x", family = "Amiri") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~target) + 
  labs(
    title = "Adumbration in Three Cases",
    subtitle = "Days since references in each case were last cited by other cases",
    x = NULL
  )

set.seed(123); bb_layout <- delete_vertices(
  graph = bb_net_positive, 
  v = which(degree(bb_net_positive) == 0)
) |> 
  #tidygraph::as_tbl_graph() |> 
  create_layout("fr")

figure5 <- bb_layout |> 
  ggraph() +
  geom_edge_fan(alpha = 2/3) +
  geom_node_point(
    mapping = aes(size = degree, fill = factor(cluster)), 
    color = "white", shape = 21, show.legend = FALSE
  ) + 
  geom_node_label(
    mapping = aes(label = str_wrap(name, 15), fill = factor(cluster), 
                  filter = name %in% principles_target), 
    size = 2, repel = TRUE, show.legend = FALSE, family = "Amiri"
  ) +
  labs(
    title = "Full Two-Mode Projection of Legal Principles", 
    subtitle = latex2exp::TeX(r"(Backbone algorithm: Stochastic Degree Sequence Model ($\alpha =$0.05))"),
    caption = "Layout algorithm: Fruchterman-Reingold"
  )


figure6 <- delete_vertices(
  graph = bb_net, 
  v = which(!V(bb_net)$name %in% case_principles_keep)
) |> 
  ggraph("graphopt")  + 
  geom_edge_fan(color = "#999999", width = 1, show.legend = FALSE) +
  geom_node_point(
    mapping = aes(fill = factor(cluster)), 
    color = "white", shape = 21, show.legend = FALSE, size = 8
  ) +
  geom_node_text(aes(label = str_wrap(name, 15)), size = 1.5) +
  labs(title = "Decision C-370-06", subtitle = "An example of principle recombination") +
  facet_edges(~ fct_rev(factor(ifelse(weight == 1, "Typical Ties", "Atypical Ties")))) 

figure7 <- stm_el |> 
  ggplot(aes(corr)) + 
  geom_histogram(color = "white") + 
  geom_rug(alpha = 1/2) + 
  labs(title = "Pairwise topic correlations", 
       subtitle = paste0(scales::comma(choose(nrow(stm_corr),2)), "coefficients"),
       x = latex2exp::TeX(r"($\rho_{ij}$)"))

set.seed(1234); corr_layout <- create_layout(corr_net, "fr")

figure8 <- corr_layout |> 
  ggraph() + 
  geom_edge_fan(aes(alpha = weight), show.legend = FALSE) +
  geom_node_point(
    mapping = aes(size = d, fill = cluster), 
    color = "white", shape = 21, show.legend = FALSE
  ) +
  geom_node_text(
    mapping = aes(label = name), 
    size = 3, show.legend = FALSE, family = "Amiri"
  ) +
  labs(title = "Positive Topic Correlations", subtitle = "Cutoff = 0.15",
       caption = "Layout algorithm: Fruchterman-Reingold")

figure9 <- case_transposition |> 
  rename(membership = weight, weight = score) |> 
  mutate(edge_type = case_when(
    sign(weight) == 1 ~ "Familiar Ties",
    sign(weight) == -1 ~ "Unfamiliar Ties"
  )) |>
  mutate(edge_type = fct_rev(edge_type)) |> 
  igraph::graph_from_data_frame(directed = FALSE) |> 
  tidygraph::as_tbl_graph() |> 
  ggraph("graphopt") + 
  geom_edge_link(
    mapping = aes(color = edge_type, label = round(weight * 100, 2), 
                  width = abs(weight)), 
    show.legend = FALSE, family = "Amiri", label_dodge = unit(20, 'mm')) +
  geom_node_point(size = 10, shape = 21, fill = "white") +
  geom_node_text(aes(label = name)) +
  facet_edges(~edge_type) +
  labs(
    caption = "Note: all edge weights have been multiplied by 100",
    subtitle = "An example of topic transposition",
    title = "Decision T-859-03"
  ) +
  scale_edge_color_manual(values = c("steelblue1", "tomato"))

appendix1 <- metadata |> 
  inner_join(cd_index) |> 
  ggplot(aes(type, cd5)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(alpha = 1/2, width = 1/4, size = 1/3) + 
  stat_summary(fun.data = mean_cl_boot, fun.args = list(conf.int = 0.9), 
               color = "red", fatten = 1) +
  labs(title = latex2exp::TeX(r"{$CD_5$ Index}"), x = NULL, y = NULL) +
  coord_cartesian(ylim = c(-1, 1))
