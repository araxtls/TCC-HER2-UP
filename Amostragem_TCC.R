# PARÂMETROS GERAIS: intervalo de confiança, variancia e poder estatístico
alpha <- 0.05
power <- 0.80

z_alpha <- qnorm(1 - alpha/2)
z_beta  <- qnorm(power)

var_kappa <- 0.05

# FUNÇÃO PARA CALCULAR O N
calc_n <- function(kappa0, kappa1, var_kappa){
  delta <- kappa1 - kappa0
  n <- ((z_alpha + z_beta)^2 * (2 * var_kappa)) / (delta^2)
  return(ceiling(n))
}

# DIVIDINDO O CALCULO DO N DE ACORDO COM OS CENÁRIOS DE TIPO TUMORAL
# 1. Alta marcação (possuem um kappa alto)
kappa0_high <- 0.75
kappa1_high <- 0.90
N_high <- calc_n(kappa0_high, kappa1_high, var_kappa)

# 2. Baixa marcação (possuem um kappa baixo)
kappa0_low <- 0.40
kappa1_low <- 0.60
N_low <- calc_n(kappa0_low, kappa1_low, var_kappa)

# RESULTADO FINAL
data.frame(
  Cenario = c("Alta concordância", "Baixa concordância"),
  kappa_semQuPath = c(kappa0_high, kappa0_low),
  kappa_comQuPath = c(kappa1_high, kappa1_low),
  N_Estimado = c(N_high, N_low)
)

N = cat("O tamanho mínimo da amostra em questão é de", sum(N_high, N_low), "lâminas"
