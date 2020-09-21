# Processo gerador de dados #

# Tamanho da amostra
t = 100

# Desvio-padrao de eta
dn = 10

# Desvio-padrao de epsilon
de = 10

# Vetor epsilon
e = rnorm(t,0,de)

# Vetor eta
n = rnorm(t,0,dn)

# Estado inicial - Parametros conhecidos
a1 = 10
P1 = 2

# Declara os vetores
alfa = rep(0,t)
y = rep(0,t)
yt = rep(0,t)

# Primeira observacao - estado inicial
alfa[1] = rnorm(1,a1,P1)

# Processo Gerador de Dados
for (i in 1:t){
  
  alfa[i+1] = alfa[i] + n[i]
  
  y[i] = alfa[i] + e[i]
  
  if (i==1){
    yt[i] = alfa[1] + e[i]
    
  }
  else{
    yt[i] = alfa[1] + sum(n[1:(i-1)]) + e[i]
    
  }
  
}

# Grafico do processo (y)
ts.plot(y)
