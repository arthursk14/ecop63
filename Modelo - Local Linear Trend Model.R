# Processo gerador de dados #
    
    # Pacote necessario para simular uma distribuicao normal multivariada
        
        require("MASS")
    
    # Tamanho da amostra
    
        t = 100
    
    # Variancia de epsilon: sigma-epsilon
        
        s.epsilon = 0.003
        H. = s.epsilon
    
    # Matriz de covariancia de eta. Como neste modelo temos dois componentes (nivel e tendencia),
    # temos duas variancias na diagonal da matriz de covariancia de eta (sigma-csi e sigma-zeta)
    
        s.csi = 0.012
        s.zeta = 1.5e-11
        
        # Coloca na matriz
            
            d = c(s.csi, s.zeta)
            
            Q. = matrix(rep(0,length(d)**2), nrow=length(d), ncol=length(d))
            
            for(i in 1:length(d)){
              for (j in 1:length(d)){
                if(i==j){
                  Q.[i,j] = d[i]
                }
              }
            }
    
    # Vetor epsilon
        
        epsilon = rnorm(t,0,H.)
    
    # Matriz eta
        
        eta = t(mvrnorm(t,rep(0,length(d)),Q.))
    
    # Estado inicial - Parametros conhecidos
    
        a1 = as.matrix(c(7.5, 0.0003))
        P1 = diag(2)*2
    
    # Declara a matriz alfa, com 2 colunas: para os componentes de nivel e tendencia.
        
        alfa = matrix(rep(0,length(d)*(t+1)), nrow=length(d), ncol=t+1)
        
    # Primeira observacao - estado inicial
        
        alfa[,1] = mvrnorm(1,a1,P1)
        
    # Cria vetor Z, que determina as variveis que impactam diretamente na equacao das observacoes, 
    # matriz R (chamada matriz de selecao, pois seleciona quais variaveis do vetor alfa sao estocasticas), e 
    # matriz T, que determina as variaveis que formam a equacao de estado, 
    # tudo conforme disposto na pagina 47 do livro
        
        Z. = c(1, 0)
        R. = matrix(c(1,0,0,1), nrow=length(d), ncol=length(d))
        T. = matrix(c(1,0,1,1), nrow=length(d), ncol=length(d))
        
    # Cria o vetor das observacoes 
        
        y = rep(0,t)
    
    # Processo Gerador de Dados
        
        for (i in 1:t){
          
          alfa[,i+1] = T. %*% alfa[,i] + R. %*% eta[,i]
          
          y[i] = Z. %*% alfa[,i] + epsilon[i]
          
        }
      
      # Grafico do processo (y)
        
          ts.plot(y)
