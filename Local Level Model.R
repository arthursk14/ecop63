### O modelo de nivel local ###
    
    # Pacote dlm, apenas para a funcao dropFirst

      require("dlm")
            
    # Dados - River Nile
            
        y = datasets::Nile
        t = length(y)
        
    # Grafico do processo (y)
        
        ts.plot(y)
        
    # Maximiza log da verossimilhanca concentrada e difusa, conforme equacao (2.63) do livro,
    # para encontrar as estimativas dos hiperparametros (sigma-epsilon, sigma-eta).
        
        logL <- function(psi){
          
          # Declara q
          q = exp(psi)
          
          # Declara os vetores
          P = rep(0,t+1)
          a = rep(0,t+1)
          F. = rep(0,t)
          v = rep(0,t)
          K = rep(0,t)
          
          # Segunda observacao
          P[2] = 1 + q
          a[2] = y[1]
          
          # Demais observacoes
          for (i in 2:t){
            
            F.[i] = P[i] + 1
            
            v[i] = y[i] - a[i]
            
            K[i] = ( P[i] / F.[i] )
            
            P[i+1] = P[i]*( 1 - K[i] ) + q
            
            a[i+1] = a[i] + K[i] * v[i]
            
          }
          
          # Estimador da variancia de epsilon, conforme equacao (2.62)
          aux = v**2 / F.
          de.hat = ( 1 / ( t - 1 ) ) * sum( aux[2:t] )
          
          # Declara a funcao de log verossimilhanca difusa e concentrada, sem as constantes do comeco (pois nao impacta a otimizacao).
          # Por padrao, a funcao optim do R minimiza o objetivo, entao os sinais estao trocados.
          
          logL = ( ( t - 1 ) / 2 ) * log(de.hat) + ( 1 / 2 ) * sum( log(F.[2:t]) )
          
          return(c(logL,de.hat))
          
        }
        
        # Otimizacao, algoritmo BFGS, conforme tabela 2.1 do livro
        
        fn <- function(p){
          return(logL(p)[1])
          
        }
        
        mlogL = optim(0, fn, method="BFGS")
        
        # Hiperparametros do modelo (sigma-epsilon e sigma-eta)
        
        de = logL( mlogL$par )[2]
        dn = de * exp( mlogL$par )
      
    # O filtro de Kalman
        
        # Inicializacao difusa, declara a1 e P1, conforme ponto 2.2.5 e secao 2.9 do livro
        
            # Valor arbitrario para a1, com a variancia tendendo à infinito, temos a2 = y1
            a1 = 0
            # Valor muito alto para variancia de a1: P1 -> Inf
            P1 = 10**7
        
        # Variancia do estimador da variavel de estado: P = Var( alfa | y )
        
            # Declara o vetor
            P = rep(0,t+1)
            
            # Primeira observacao - Estado inicial
            P[1] = P1
            
            # Demais observacoes
            for (i in 1:t){
              P[i+1] = ( ( P[i] * de ) / ( P[i] + de ) ) + dn
              
            }
        
        # Estimador da variavel de estado: a = E( alfa | y )
            
            # Declara o vetor
            a = rep(0,t+1)
            
            # Primeira observacao - Estado inicial
            a[1] = a1
            
            # Demais observacoes
            for (i in 1:t){
              a[i+1] = a[i] + ( P[i] / (P[i] + de) ) * ( ( y[i] - a[i] )  )
              
            }
            
        # Coloca a serie do modelo estimado no grafico
        ts.plot(cbind(y,dropFirst(a)),col=c("black","red"))
      
        # Monta serie do erro de previsao (v), da variancia do erro de previsao (F) 
        # e do ganho de Kalman (K). Todos ja estavam implicitos nas recursoes acima.
        
            # Declara os vetores
            v = rep(0,t)
            F. = rep(0,t)
            K = rep(0,t)
            
            # Coloca as observacoes
            for (i in 1:t){
              v[i] = y[i] - a[i]
              
              F.[i] = P[i] + de
              
              K[i] = P[i] / F.[i]
              
            }
            
    # Suavizacao, recursao reversa, para estimar a variavel de estado com base em toda a informacao realizada,
    # conforme secao 2.4 do livro
    
        # Notacoes auxiliares: L, conforme equacao 2.30; r, conforme equacao 2.36; e N, conforme equacao 2.42
          
            # Declara os vetores
            L = rep(0,t)
            r = rep(0,t)
            N = rep(0,t)
            
            # Coloca as observacoes
            for (i in t:1){
              L[i] = ( 1 - K[i] )
              
              if(i==1){
                # Primeira observacao: Pela construcao da recursao reversa, fica fora do vetor
                r_0 = ( v[i] / F.[i] ) + ( L[i] * r[i] )
                N_0 = ( 1 / F.[i] ) + ( ( L[i]**2 ) * N[i] ) 
                
              }
              else{
                # Tambem pela construcao, a ultima observacao fica zerada
                r[i-1] = ( v[i] / F.[i] ) + ( L[i] * r[i] )
                N[i-1] = ( 1 / F.[i] ) + ( ( L[i]**2 ) * N[i] )
                
              }
              
            }
            
        # Calcula a variavel de estado suavizada, conforme equacao 2.37
        
            # Declara o vetor alfa.hat
            alfa.hat = rep(0,t)
            
            # Coloca as observacoes
            for (i in t:1){
              if(i==1){
                alfa.hat[i] = a[i] + ( P[i] * r_0 )
                
              }
              else{
                alfa.hat[i] = a[i] + ( P[i] * r[i-1] )
                
              }
              
            }
            
        # Calcula a variancia da variavel de estado suavizada, conforme equacao 2.43
            
            # Declara o vetor V
            V = rep(0,t)
            
            # Coloca as observacoes
            for (i in t:1){
              if(i==1){
                v[i] = P[i] - ( ( P[i]**2 ) * N_0 )
                
              }
              else{
                v[i] = P[i] - ( ( P[i]**2 ) * N[i-1] )
                
              }
              
            } 
                
    # Coloca a serie do modelo estimado (suavizado) no grafico
    ts.plot(cbind(y,(dropFirst(a)),(alfa.hat)),col=c("black","red","blue"))
