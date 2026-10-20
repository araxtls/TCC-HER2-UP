## PACOTES NECESSÁRIOS
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(irr)
library(psych)
library(krippendorff) ##rever
library(writexl)
library(ggplot2)

dados <- read_csv(
  "/Users/rafaellatoledoaranha/Desktop/QUPATH_HER2.csv",
  show_col_types = FALSE
)
glimpse(dados)

categorias <- c("0", "0+", "1+", "2+") ## transformando categorias em valores ordinais
categorias2 <- c("0", "0+", "1+")

comparacao <- data.frame( ## data frame com os dados mais importantes renomeados
  ID = dados[["ID IMUNO"]],
  
  APMS_MO = trimws(as.character(
    dados[["APMS \n(SEM QUPATH)"]]
  )),
  
  APP_MO = trimws(as.character(
    dados[["APP \n(SEM QUPATH)"]]
  )),
  
  APMS_QuPath = trimws(as.character(
    dados[["APMS \n(COM QUPATH)"]]
  )),
  
  APP_QuPath = trimws(as.character(
    dados[["APP \n(COM QUPATH)"]]
  )),
  
  DX_QuPath = trimws(as.character(
    dados[["Qupath Dx"]]
  ))
)

comparacao[comparacao == ""] <- NA ## descartando células vazias (não é para ter nenhuma)
lapply(comparacao[-1], table, useNA = "ifany") ## conferindo categorias encontradas

calcular_concordancia <- function(avaliadora1, avaliadora2, nome) { ## funcao para calcular os kappas
  
    completos <- complete.cases(avaliadora1, avaliadora2) ## mantendo apenas casos completos
  
  a <- avaliadora1[completos]
  b <- avaliadora2[completos]
  
  a <- factor(a, levels = categorias, ordered = TRUE)
  b <- factor(b, levels = categorias, ordered = TRUE)
  
  matriz <- table(   ## matriz de confusão
    Avaliadora1 = a,
    Avaliadora2 = b
  )
  
  concordancia <- mean(a == b) * 100 ## concordancia exata
  
  kappa_simples <- irr::kappa2( ## kappa de cohen
    data.frame(a, b),
    weight = "unweighted"
  )
  
  kappa_linear <- irr::kappa2( ## kappa de cohen ponderado linearmente
    data.frame(a, b),
    weight = "equal"
  )
  
  kappa_quadratico <- irr::kappa2( ## kappa de cohen ponderado quadraticamente
    data.frame(a, b),
    weight = "squared"
  )
  
  cat("\n========================================\n") ## resultados
  cat(nome, "\n")
  cat("========================================\n")
  
  cat("Casos totais:", length(avaliadora1), "\n")
  cat("Casos completos:", sum(completos), "\n")
  cat("Concordância exata:", round(concordancia, 2), "%\n")
  
  cat("\nMatriz de confusão:\n")
  print(matriz)
  
  cat("\nKappa não ponderado:\n")
  print(kappa_simples)
  
  cat("\nKappa ponderado linearmente:\n")
  print(kappa_linear)
  
  cat("\nKappa ponderado quadraticamente:\n")
  print(kappa_quadratico)
  
  return(list(
    nome = nome,
    n = sum(completos),
    concordancia_exata = concordancia,
    matriz = matriz,
    kappa_simples = kappa_simples,
    kappa_linear = kappa_linear,
    kappa_quadratico = kappa_quadratico
  ))
}

resultado_MO <- calcular_concordancia( ## concordancia entre as médicas - momento 1
  comparacao$APMS_MO,
  comparacao$APP_MO,
  "Concordância entre avaliadoras — microscopia óptica"
)

resultado_QuPath <- calcular_concordancia(  ## concordancia entre as médicas - momento 2
  comparacao$APMS_QuPath,
  comparacao$APP_QuPath,
  "Concordância entre avaliadoras — com QuPath"
)

resultado_intraobs_APMS <- calcular_concordancia( ## concordancia intraobservador APMS
  comparacao$APMS_MO,
  comparacao$APMS_QuPath,
  "Concordância intra avaliadora APMS"
)

resultado_intraobs_APP <- calcular_concordancia( ## concordancia intraobservador APP
  comparacao$APP_MO,
  comparacao$APP_QuPath,
  "Concordância intra avaliadora APP"
)

validade_qupath_sozinho <- calcular_concordancia( ## concordancia algoritmo - humano
  comparacao$APMS_QuPath,
  comparacao$DX_QuPath,
  "Concordancia intra avaliadora APMS - QuPath"
)

validade_qupath_sozinho <- calcular_concordancia( ## concordancia algoritmo - humano
  comparacao$APP_QuPath,
  comparacao$DX_QuPath,
  "Concordancia intra avaliadora APP - QuPath"
)

reclassificacoes_APMS <- ( ## taxa de reclassificacao APMS
  comparacao$APMS_MO[validos_APMS] !=
    comparacao$APMS_QuPath[validos_APMS]
)

n_reclassificados_APMS <- sum(reclassificacoes_APMS)
n_validos_APMS <- sum(validos_APMS)

pct_reclassificados_APMS <-
  100 * n_reclassificados_APMS / n_validos_APMS

reclassificacoes_APP <- (  ## taxa de reclassificacao APMS
  comparacao$APP_MO[validos_APP] !=
    comparacao$APP_QuPath[validos_APP]
)

n_reclassificados_APP <- sum(reclassificacoes_APP)
n_validos_APP <- sum(validos_APP)

pct_reclassificados_APP <-
  100 * n_reclassificados_APP / n_validos_APP

resumo_reclassificacao <- data.frame( ## resumo
  Avaliadora = c("APMS", "APP"),
  Casos_avaliados = c(n_validos_APMS, n_validos_APP),
  Casos_reclassificados = c(
    n_reclassificados_APMS,
    n_reclassificados_APP
  ),
  Percentual_reclassificado = round(
    c(pct_reclassificados_APMS, pct_reclassificados_APP),
    2
  )
)
print(resumo_reclassificacao)

transicao_APMS <- table( ## matriz de transição APMS
  MO = factor(
    comparacao$APMS_MO,
    levels = categorias,
    ordered = TRUE
  ),
  QuPath = factor(
    comparacao$APMS_QuPath,
    levels = categorias,
    ordered = TRUE
  )
)

transicao_APP <- table( ## matriz de transição APP
  MO = factor(
    comparacao$APP_MO,
    levels = categorias,
    ordered = TRUE
  ),
  QuPath = factor(
    comparacao$APP_QuPath,
    levels = categorias,
    ordered = TRUE
  )
)

cat("\nTransição de escores — APMS:\n")
print(transicao_APMS)

cat("\nTransição de escores — APP:\n")
print(transicao_APP)

## MAPAS DE CALOR
tema_heatmap <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "plain",
      size = 14,
      hjust = 0.5
    ),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.title = element_text(face = "plain"),
    legend.text = element_text(face = "plain")
  )

matriz_MO <- table( ## APMS - APP (microscopia optica)
  factor(comparacao$APMS_MO, levels = categorias2),
  factor(comparacao$APP_MO, levels = categorias2)
)
dados_heatmap_MO <- as.data.frame(matriz_MO)

colnames(dados_heatmap_MO) <- c(
  "APMS", "APP", "Frequencia"
)

grafico_heatmap_MO <- ggplot(
  dados_heatmap_MO,
  aes(x = APP, y = APMS, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias2,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias2),
    drop = FALSE
  ) +
  labs(
    title = "Concordância entre observadoras — MO",
    x = "Observador A",
    y = "Observador B"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_MO)

matriz_QuPath <- table( ## APMS - APP (analise digital)
  factor(comparacao$APMS_QuPath, levels = categorias),
  factor(comparacao$APP_QuPath, levels = categorias)
)
dados_heatmap_QuPath <- as.data.frame(matriz_QuPath)

colnames(dados_heatmap_QuPath) <- c(
  "APMS", "APP", "Frequencia"
)

grafico_heatmap_QuPath <- ggplot(
  dados_heatmap_QuPath,
  aes(x = APP, y = APMS, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Concordância entre observadoras — QuPath",
    x = "Observador A",
    y = "Observador B"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_QuPath)

matriz_APMS <- table(  ## APMS - antes e após qupath
  factor(comparacao$APMS_MO, levels = categorias),
  factor(comparacao$APMS_QuPath, levels = categorias)
)

dados_heatmap_APMS <- as.data.frame(matriz_APMS)

colnames(dados_heatmap_APMS) <- c(
  "MO", "QuPath", "Frequencia"
)

grafico_heatmap_APMS <- ggplot(
  dados_heatmap_APMS,
  aes(x = QuPath, y = MO, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Reclassificação do Observador A",
    x = "Análise digital",
    y = "Microscopia óptica"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_APMS)

matriz_APP <- table(   ## APP - antes e após qupath
  factor(comparacao$APP_MO, levels = categorias),
  factor(comparacao$APP_QuPath, levels = categorias)
)
dados_heatmap_APP <- as.data.frame(matriz_APP)

colnames(dados_heatmap_APP) <- c(
  "MO", "QuPath", "Frequencia"
)

grafico_heatmap_APP <- ggplot(
  dados_heatmap_APP,
  aes(x = QuPath, y = MO, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Reclassificação do Observadora B",
    x = "Análise digital",
    y = "Microscopia óptica"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_APP)

length(comparacao$APP_MO)
length(comparacao$DX_QuPath)
nrow(comparacao)

matriz_QuPath0 <- table(   ## APP MO x QuPath
  factor(comparacao$APP_MO, levels = categorias),
  factor(comparacao$DX_QuPath, levels = categorias)
)
dados_heatmap_QuPath0 <- as.data.frame(matriz_QuPath0)

colnames(dados_heatmap_QuPath0) <- c(
  "MO APP", "QuPath", "Frequencia"
)

grafico_heatmap_QuPath0 <- ggplot(
  dados_heatmap_QuPath0,
  aes(x = QuPath, y = `MO APP`, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Concordância APP (MO) - QuPath",
    x = "Análise digital isolada",
    y = "Observador B (MO)"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_QuPath0)

matriz_QuPath1 <- table(   ## APP QuPath x QuPath
  factor(comparacao$APP_QuPath, levels = categorias),
  factor(comparacao$DX_QuPath, levels = categorias)
)
dados_heatmap_QuPath1 <- as.data.frame(matriz_QuPath1)

colnames(dados_heatmap_QuPath1) <- c(
  "QuPath APP", "QuPath", "Frequencia"
)

grafico_heatmap_QuPath1 <- ggplot(
  dados_heatmap_QuPath1,
  aes(x = QuPath, y = `QuPath APP`, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Concordância APP (QuPath) - QuPath",
    x = "Análise digital isolada",
    y = "Observador B (Análise digital)"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_QuPath1)

matriz_QuPath2 <- table(   ## APMS QuPath x QuPath
  factor(comparacao$APMS_QuPath, levels = categorias),
  factor(comparacao$DX_QuPath, levels = categorias)
)
dados_heatmap_QuPath2 <- as.data.frame(matriz_QuPath2)

colnames(dados_heatmap_QuPath2) <- c(
  "QuPath APMS", "QuPath", "Frequencia"
)

grafico_heatmap_QuPath2 <- ggplot(
  dados_heatmap_QuPath2,
  aes(x = QuPath, y = `QuPath APMS`, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Concordância APMS (QuPath) - QuPath",
    x = "Análise digital isolada",
    y = "Observador A (Análise digital)"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_QuPath2)

matriz_QuPath3 <- table(   ## APMS MO x QuPath
  factor(comparacao$APMS_MO, levels = categorias),
  factor(comparacao$DX_QuPath, levels = categorias)
)
dados_heatmap_QuPath3 <- as.data.frame(matriz_QuPath3)

colnames(dados_heatmap_QuPath3) <- c(
  "MO APMS", "QuPath", "Frequencia"
)

grafico_heatmap_QuPath3 <- ggplot(
  dados_heatmap_QuPath3,
  aes(x = QuPath, y = `MO APMS`, fill = Frequencia)
) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(
    aes(label = Frequencia),
    size = 5,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#FFF5F8",
    high = "#FF82AB",
    name = "Número de casos"
  ) +
  scale_x_discrete(
    limits = categorias,
    drop = FALSE
  ) +
  scale_y_discrete(
    limits = rev(categorias),
    drop = FALSE
  ) +
  labs(
    title = "Concordância APMS (MO) - QuPath",
    x = "Análise digital isolada",
    y = "Observador A (MO)"
  ) +
  coord_fixed() +
  tema_heatmap
print(grafico_heatmap_QuPath3)

## GRAFICOS DE BARRAS AGRUPADAS
dados_comparacao <- comparacao %>%
  select(APMS_MO, APP_MO, APMS_QuPath, APP_QuPath) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variavel",
    values_to = "Escore"
  ) %>%
  mutate(
    Observadora = case_when(
      grepl("^APMS", Variavel) ~ "Observador A",
      grepl("^APP", Variavel) ~ "Observador B"
    ),
    Momento = case_when(
      grepl("_MO$", Variavel) ~ "Microscopia óptica",
      grepl("_QuPath$", Variavel) ~ "Análise digital"
    ),
    Escore = as.character(Escore),
    Escore = factor(
      Escore,
      levels = categorias
    )
  )

# Frequências por método, respeitando as categorias de cada um
frequencias_comparacao <- dados_comparacao %>%
  count(Momento, Observadora, Escore, .drop = FALSE) %>%
  filter(
    (Momento == "Microscopia óptica" &
       Escore %in% categorias2) |
      Momento == "Análise digital"
  ) %>%
  mutate(
    Escore = factor(
      as.character(Escore),
      levels = categorias
    ),
    Momento = factor(
      Momento,
      levels = c("Microscopia óptica", "Análise digital")
    )
  )

cores_observadoras <- c(
  "Observador A" = "#FFB5C5",
  "Observador B" = "#FF82AB"
)

grafico_comparativo <- ggplot(
  frequencias_comparacao,
  aes(x = Escore, y = n, fill = Observadora)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(
    ~ Momento,
    nrow = 1,
    scales = "free_x"
  ) +
  scale_fill_manual(
    values = cores_observadoras,
    limits = c("Observador A", "Observador B"),
    breaks = c("Observador A", "Observador B"),
    drop = FALSE
  ) +
  scale_x_discrete(
    limits = function(x) {
      x[x %in% categorias]
    },
    drop = TRUE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Distribuição dos Escores de HER2 por Observador e Método",
    x = "Escore HER2",
    y = "Número de casos",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "plain",
      size = 14,
      hjust = 0.5
    ),
    axis.title = element_text(face = "plain"),
    strip.text = element_text(
      face = "plain",
      size = 12
    ),
    legend.title = element_text(face = "plain"),
    legend.text = element_text(face = "plain"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
print(grafico_comparativo)
