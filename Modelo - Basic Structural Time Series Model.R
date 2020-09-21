# Processo gerador de dados #
    
    # Pacote necessario para simular uma distribuicao normal multivariada
        
        require("MASS")
    
    # Tamanho da amostra
    
        t = 100
    
    # Variancia de epsilon: sigma-epsilon
        
        H = 1
    
    # Matriz de covariancia de eta. Como neste modelo temos tres componentes (nivel, tendencia e sazonalidade),
    # temos tres variancias na diagonal da matriz de covariancia eta (sigma-csi, sigma-zeta e sigma-omega),
    # conforme disposto na pagina 47 do livro
    
        s.csi = 1
        s.zeta = 1
        s.omega = 1
        
        # Coloca na matriz
            
            d = c(s.csi, s.zeta, s.omega)
            
            Q = matrix(rep(0,3**2), nrow=3, ncol=3)
            
            for(i in 1:3){
              for (j in 1:3){
                if(i==j){
                  Q[i,j] = d[i]
                }
              }
            }
    
    # Vetor epsilon
        
        epsilon = rnorm(t,0,H)
    
    # Matriz eta
        
        eta = t(mvrnorm(t,c(0,0,0),Q))
    
    # Estado inicial - Parametros conhecidos
    
        a1 = t(c(0, 0, 0, 0, 0))
        
        d = c(1, 1, 1, 1, 1)
        P1 = matrix(rep(0,5**2), nrow=5, ncol=5)
        
        for(i in 1:5){
          for (j in 1:5){
            if(i==j){
              P1[i,j] = d[i]
            }
          }
        }
    
    # Declara a matriz alfa, com 5 colunas: 2 para os componentes de nivel e tendencia.
    # 3 para o compontente sazonal: (s-1) = 3 elementos, conforme disposto na pagina 47 do livro.
        
        alfa = matrix(rep(0,5*(t+1)), nrow=5, ncol=t+1)
        
    # Primeira observacao - estado inicial
        
        alfa[,1] = mvrnorm(1,a1,P1)
        
    # Cria vetor Z, que determina as variveis que impactam diretamente na equacao das observacoes, 
    # matriz R (chamada matriz de selecao, pois seleciona quais variaveis do vetor alfa sao estocasticas), e 
    # matriz T, que determina as variaveis que formam a equacao de estado, 
    # tudo conforme disposto na pagina 47 do livro
        
        Z = c(1, 0, 1, 0, 0)
        R = matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0), nrow=5, ncol=3)
        T. = matrix(c(1,0,0,0,0,1,1,0,0,0,0,0,-1,1,0,0,0,-1,0,1,0,0,-1,0,0), nrow=5, ncol=5)
        
    # Cria o vetor das observacoes 
        
        y = rep(0,t)
    
    # Processo Gerador de Dados
        
        for (i in 1:t){
          
          alfa[,i+1] = T.%*%alfa[,i] + R%*%eta[,i]
          
          y[i] = Z%*%alfa[,i] + epsilon[i]
          
        }
      
      # Grafico do processo (y)
        
          ts.plot(y)
