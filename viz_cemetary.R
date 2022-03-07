## PROTOTYPES

## RUN THE FIRST FEW CHUNKS OF THE second-year-paper.Rmd file first

# ADUMBRATION -------------------------------------------------------
metadata |> 
  ggplot(aes(adum_max / 365, y = factor(year))) + 
  geom_boxplot(outlier.alpha = 1/4, outlier.size = 1/2) +
  stat_summary(fun.data = mean_cl_boot, fatten = 1, color = "red") +
  labs(y = NULL, title = "Adumbration by year", x = "years",
       subtitle = "Largest amount of time since last cited")

adumbration |> 
  mutate(year = factor(ccc:::extract_year(target))) |> 
  ggplot(aes(as.numeric(days / 365), y = year)) + 
  geom_boxplot(outlier.alpha = 1/4, outlier.size = 1/2) +
  stat_summary(fun.data = mean_se, fatten = 1, color = "red") +
  labs(y = NULL, title = "Adumbration by year", x = "years",
       subtitle = "Largest amount of time since last cited")

quibble <- function(x, q = c(0.9, 0.95, 0.99)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}

adumbration |> 
  mutate(year = factor(ccc:::extract_year(target))) |> 
  group_by(year) |> 
  summarize(
    quibble(days)    
  ) |> 
  mutate(days_q = factor(days_q)) |> 
  ggplot(aes(as.numeric(days / 365), year, color = days_q)) + 
  geom_point()

# RECOMBINATION -------------------------------------------------------

metadata |> 
  pivot_longer(c(atypical, typical), names_to = "recombination", values_to = "score") |> 
  ggplot(aes(score, recombination)) +
  geom_jitter(height = 1/8, alpha = 1/4, size = 1/5) +
  geom_violin(alpha = 1/2, show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, fatten = 1, show.legend = FALSE, color = "red") +
  scale_x_log10()

progresividad_target <- c("principio de progresividad", "principio pro homine", "principio in dubio pro operario", "principio de universalidad", "principio de eficiencia del tributo", "principio de proporcionalidad", "principio de certeza")

igraph::make_ego_graph(
  graph = delete_vertices(
    graph = bb_net_positive, 
    v = which(degree(bb_net_positive) == 0)
  ),
  nodes = "principio de progresividad"
) |> 
  pluck(1) |> 
  ggraph("kk") + 
  geom_edge_fan(alpha = 2/3) +
  geom_node_point(aes(size = degree, fill = factor(cluster)), 
                  color = "white", shape = 21, show.legend = FALSE) +
  geom_node_label(
    mapping = aes(label = str_wrap(name, 20), fill = factor(cluster),
                  filter = name %in% progresividad_target),
    size = 2, family = "Amiri", show.legend = FALSE) +
  labs(title = "Principle of Progressivity", subtitle = "Ego Network",
       caption = "Layout algorithm: Kamada-Kawai") 

igraph::make_ego_graph(
  graph = delete_vertices(
    graph = bb_net, 
    v = which(degree(bb_net) == 0)
  ),
  nodes = "principio de progresividad"
) |> 
  pluck(1) |> 
  ggraph("graphopt")  + 
  geom_edge_fan(alpha = 2/3) +
  geom_node_point(
    mapping = aes(size = degree, fill = factor(cluster)), 
    color = "white", shape = 21, show.legend = FALSE
  ) +
  labs(title = "Principle of Progressivity", subtitle = "Ego Network",
       caption = "Layout algorithm: Kamada-Kawai") +
  facet_edges(~factor(sign))


## Transposition
metadata |> 
  pivot_longer(c(familiar, unfamiliar), names_to = "transposition", values_to = "score") |> 
  ggplot(aes(score, y = transposition, group = transposition)) + 
  geom_jitter(height = 1/8, alpha = 1/4, size = 1/5) +
  geom_violin(alpha = 1/2, show.legend = FALSE) + 
  stat_summary(fun.data = mean_se, fatten = 1, show.legend = FALSE, color = "red") +
  labs(y = NULL, title = "Transposition Raw Scores") 

## DESCRIPTIVE

tidyr::crossing(
  x = seq(0, 5, 1),
  year = seq(1992, 2021)
) |> 
  mutate(prop = map2_dbl(x, year, function(x, y) {
    mean(metadata$indegree[which(metadata$year == y)] <= x)
  })) |> 
  group_by(year) |> 
  mutate(l = scales::percent(max(prop))) |> 
  ggplot(aes(x, prop)) + 
  geom_line() + 
  geom_text(aes(label = l), x = 2.5, y = 0.5) +
  facet_wrap(~factor(year)) 
scale_x_continuous(labels = seq(0, 5, 1), breaks = seq(0, 5, 1)) +
  labs(x = "indegree", title = "Cumulative Distribution of All Citations by Year")