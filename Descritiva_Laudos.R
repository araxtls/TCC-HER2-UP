library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)

dados <- read_csv( ## acessando a tabela
  "/Users/rafaellatoledoaranha/Desktop/Estatística Descritiva TCC.csv")

glimpse(dados)
dim(dados)
names(dados)

str(dados[c("ER%", "PR%", "KI-67")])
head(dados[c("ER%", "PR%", "KI-67")], 10)

dados <- dados %>% ## arrumando os dados representados com %
  mutate(
    across(
      c(`ER%`, `PR%`, `KI-67`),
      ~ .x / 100
    )
  )

dados <- dados %>%
  mutate(
    across(
      c(`ER%`, `PR%`, `KI-67`),
      ~ parse_number(as.character(.x))
    )
  )

dados <- dados %>% ## categorias
  mutate(
    HER2 = factor(
      `CLASSIFICAÇÃO HER2`,
      levels = c("0", "0+", "1+"),
      ordered = TRUE
    ),
    
    STATUS_HER2 = factor(
      STATUS_HER2,
      levels = c("Negative", "Ultralow", "Low")
    ),
    
    PERFIL_MOLECULAR = factor(
      PERFIL_MOLECULAR,
      levels = c(
        "Luminal A",
        "Luminal B",
        "Triplo negativo"
      )
    ),
    
    ER_percent = `ER%` * 100,
    PR_percent = `PR%` * 100,
    KI67_percent = `KI-67` * 100
  )
summary(dados[c("ER%", "PR%", "KI-67")])
summary(dados[c("ER_percent", "PR_percent", "KI67_percent")])

variaveis_categoricas <- c(
  "DIAGNÓSTICO_HISTOPATOLÓGICO",
  "CARCINOMA_IN_SITU",
  "NECROSE_TUMORAL",
  "GRAU_IN_SITU",
  "MICROCALCIFICAÇÕES",
  "SUBTIPO_HISTOLÓGICO_IN_SITU",
  "PLEOMORFISMO_NUCLEAR",
  "FORMAÇÃO_TUBULAR",
  "CONTAGEM_MITÓTICA",
  "GRAU_HISTOLÓGICO",
  "INVASÃO_ANGIOLINFÁTICA",
  "INVASÃO PERINEURAL",
  "HER2",
  "STATUS_HER2",
  "PERFIL_MOLECULAR"
)

tabela_categoricas <- map_dfr( ## frequencias das variaveis categoricas
  variaveis_categoricas,
  function(var) {
    
    dados %>%
      count(.data[[var]], name = "n") %>%
      mutate(
        variavel = var,
        percentual = round(100 * n / sum(n), 1)
      ) %>%
      select(
        variavel,
        categoria = all_of(var),
        n,
        percentual
      )
  }
)
tabela_categoricas
print(tabela_categoricas, n = Inf)

dados <- dados %>% ## resumo variaveis numericas
  mutate(
    infiltrado_percent = case_when(
      `INFILTRADO INFLAMATORIO ESTROMAL ESTIMADO` == "<5" ~ 5,
      TRUE ~ as.numeric(
        `INFILTRADO INFLAMATORIO ESTROMAL ESTIMADO`
      )
    )
  )

resumo_quantitativas <- dados %>%
  summarise(
    
    # ER
    ER_media = mean(ER_percent, na.rm = TRUE),
    ER_DP = sd(ER_percent, na.rm = TRUE),
    ER_mediana = median(ER_percent, na.rm = TRUE),
    ER_Q1 = quantile(ER_percent, 0.25, na.rm = TRUE),
    ER_Q3 = quantile(ER_percent, 0.75, na.rm = TRUE),
    
    # PR
    PR_media = mean(PR_percent, na.rm = TRUE),
    PR_DP = sd(PR_percent, na.rm = TRUE),
    PR_mediana = median(PR_percent, na.rm = TRUE),
    PR_Q1 = quantile(PR_percent, 0.25, na.rm = TRUE),
    PR_Q3 = quantile(PR_percent, 0.75, na.rm = TRUE),
    
    # Ki-67
    KI67_media = mean(KI67_percent, na.rm = TRUE),
    KI67_DP = sd(KI67_percent, na.rm = TRUE),
    KI67_mediana = median(KI67_percent, na.rm = TRUE),
    KI67_Q1 = quantile(KI67_percent, 0.25, na.rm = TRUE),
    KI67_Q3 = quantile(KI67_percent, 0.75, na.rm = TRUE),
    
    # Infiltrado inflamatório
    Infiltrado_media = mean(infiltrado_percent, na.rm = TRUE),
    Infiltrado_DP = sd(infiltrado_percent, na.rm = TRUE),
    Infiltrado_mediana = median(infiltrado_percent, na.rm = TRUE),
    Infiltrado_Q1 = quantile(infiltrado_percent, 0.25, na.rm = TRUE),
    Infiltrado_Q3 = quantile(infiltrado_percent, 0.75, na.rm = TRUE)
  )
print(resumo_quantitativas, width = Inf)

dados_boxplot <- dados %>% ## boxplot para variaveis numericas (marcadores que não o HER2)
  select(ER_percent, PR_percent, KI67_percent) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variavel",
    values_to = "valor"
  ) %>%
  mutate(
    variavel = recode(
      variavel,
      ER_percent = "ER (%)",
      PR_percent = "PR (%)",
      KI67_percent = "Ki-67 (%)"
    )
  )

ggplot(dados_boxplot, aes(x = variavel, y = valor)) +
  geom_boxplot(
    fill = "#FFB5C5",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.5
  ) +
  scale_x_discrete(
    limits = c("ER (%)", "PR (%)", "Ki-67 (%)")
  ) +
  labs(
    x = NULL,
    y = "Expressão (%)"
  ) +
  theme_classic()

grafico_barras <- function(data, variavel, titulo, eixo_x = NULL) { ## grafico de barras para os principais CA invasivos
  
  tabela <- data %>%
    count({{ variavel }}, name = "n") %>%
    mutate(
      percentual = 100 * n / sum(n),
      rotulo = paste0(n, " (", round(percentual, 1), "%)")
    )
  
  ggplot(tabela, aes(x = {{ variavel }}, y = n)) +
    geom_col(fill = "#FFB5C5", width = 0.7) +
    geom_text(
      aes(label = rotulo),
      vjust = -0.4,
      size = 4
    ) +
    labs(
      title = titulo,
      x = eixo_x,
      y = "Número de casos"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5
      ),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    expand_limits(y = max(tabela$n) * 1.15)
}

dados_graficos <- dados %>%
  mutate(
    HER2 = factor(
      HER2,
      levels = c("0", "0+", "1+"),
      ordered = TRUE
    ),
    
    DIAGNÓSTICO_HISTOPATOLÓGICO = factor(
      DIAGNÓSTICO_HISTOPATOLÓGICO,
      levels = c(
        "Carcinoma ductal invasor",
        "Carcinoma lobular invasor",
        "Carcinoma tubular invasor",
        "Carcinoma túbulo - lobular invasor",
        "Carcinoma micropapilar invasor"
      ),
      labels = c(
        "Ductal sem tipo
específico",
        "Lobular",
        "Tubular",
        "Túbulo-lobular",
        "Micropapilar"
      )
    )
  )

grafico_barras(
  dados_graficos,
  DIAGNÓSTICO_HISTOPATOLÓGICO,
  "SUBTIPOS HISTOLÓGICOS",
  "Diagnóstico"
)

dados_cis <- data.frame( ## grafico barras ausencia / presenca de CA in situ
  categoria = c(
    "Ausência",
    "Presença",
    "Presença",
    "Presença"
  ),
  subtipo = c(
    "Ausência",
    "Cribriforme",
    "Sólido",
    "Comedocarcinoma"
  ),
  n = c(42, 10, 5, 3)
) %>%
  mutate(
    categoria = factor(
      categoria,
      levels = c("Ausência", "Presença")
    ),
    subtipo = factor(
      subtipo,
      levels = c(
        "Ausência",
        "Cribriforme",
        "Sólido",
        "Comedocarcinoma"
      )
    )
  )

dados_cis <- dados_cis %>%
  mutate(
    percentual_amostra = 100 * n / 60,
    percentual_presenca = ifelse(
      categoria == "Presença",
      100 * n / 18,
      NA_real_
    ),
    rotulo = ifelse(
      categoria == "Ausência",
      paste0(n, "\n(70,0%)"),
      paste0(
        n, "\n(",
        sprintf("%.1f", percentual_presenca),
        "%)"
      )
    )
  )

ggplot( ## grafico de barras empilhadas para subtipos de CA in situ
  dados_cis,
  aes(x = categoria, y = n, fill = subtipo)
) +
  geom_col(
    width = 0.65,
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    aes(label = rotulo),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Ausência" = "#FF82AB",
      "Cribriforme" = "#FFB5C5",
      "Sólido" = "#EEA9B8",
      "Comedocarcinoma" = "#CD919E"
    ),
    breaks = c(
      "Ausência",
      "Cribriforme",
      "Sólido",
      "Comedocarcinoma"
    ),
    labels = c(
      "Ausência",
      "Cribriforme",
      "Sólido",
      "Comedocarcinoma"
    ),
    name = "Categoria histológica"
  ) +
  scale_x_discrete(
    labels = c(
      "Ausência",
      "Presença"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0, 60, 10),
    limits = c(0, 60),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Carcinoma in situ na amostra",
    x = NULL,
    y = "Número de casos"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "plain",
      size = 14,
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title = element_text(face = "plain"),
    legend.text = element_text(size = 10)
  )
