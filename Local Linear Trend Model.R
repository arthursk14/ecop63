### O modelo estrutural basico para series temporais ###
    
    # Pacote dlm, para a funcao dropFirst e para estimar os hiperparametros por maxima verossimilhanca
  
        require("dlm")
        
    # Pacote quantmod, para pegar os dados dos preços das ações do Yahoo Finance
        
        require("quantmod")
        
    # Dados - UKDriverDeaths
        
        y = log(datasets::UKDriverDeaths)
        t = length(y)
        
    # Dados - Alphabet
        
        # Declara as datas
        
            d.ini = as.Date("2007-01-03")
            d.fim = as.Date("2020-09-17")
        
        # Pega dados
            
            getSymbols("GOOG", from=d.ini, to=d.fim)
        
        # Monta a série de preços de fechamento - ultimos 100 dias
            
            final = length(GOOG[,1])
            y = GOOG[(final-100):final,4]
      
    # Grafico do processo (y)
            
        t = length(y)
    
        ts.plot(y)
        
    # Matrizes que determinam o modelo de nivel local com tendencia estocastica
        
        Z. = matrix(c(1, 0), nrow = 1, ncol = 2)
        T. = matrix(c(1,0,1,1), nrow = 2, ncol = 2)
        R. = diag(2)
        
    # Hiperparametros do modelo (sigma-epsilon, sigma-csi e sigma-zeta)
        
        # Estima por maxima verossimilhanca
        
            # Modelo
        
                ordem = 2 # Representa um modelo de nivel com tendencia estocastica
                fn <- function(p){
                  dlmModPoly(order = ordem, dV = exp(p[1]), dW = exp(p[2:3]),
                             m0 = c(0, 0), C0 = 1e7*diag(2))
                }
                
            # Encontra os hiperparametros que maximizam a funcao de log verossimilhanca para o modelo
                
                fit = dlmMLE(y, rep(1,3), fn)
                mod = fn(fit$par)
                
            # Hiperparametros - Variancias - Observacao (Componente Irregular) e Estado
                
                s.epsilon = V(mod)
                s.csi = W(mod)[1]
                s.zeta = W(mod)[2]
        
        # Matrizes de covariancia
        
            H. = s.epsilon
            Q. = matrix(c(s.csi,0,0,s.zeta), nrow = 2, ncol = 2)
        
    # O filtro de Kalman
            
        # Declara o vetor y.hat e o vetor v, dos erros de previsao
            
            y.hat = rep(0,t+1)
            v = rep(0,t)
            
        # Declara os arrays para salvar as matrizes que dependem de t
            
            a = array(dim = c(2,1,t+1))
            P = array(dim = c(2,2,t+1))
            
            F. = rep(0,t)
            K. = array(dim = c(2,1,t))
            L. = array(dim = c(2,2,t))
        
        # Inicializacao difusa aproximada, conforme equacao (5.4) do livro
            
            a[,,1] = matrix(c(0, 0), nrow = 2, ncol = 1)
            P[,,1] = diag(2)*1e7 + R. %*% Q. %*% t(R.) 
            
        
        # Recursao para estimar y para i = 1, ..., t por meio do filtro, conforme equacoes (4.24)
            
            y.hat[1] = Z. %*% matrix(a[,,1])
            
            for (i in 1:t){
              
              # Erro de previsao um passo a frente
              
                  v[i] = y[i] - Z. %*% matrix(a[,,i])
              
              # Notacao adicional: F, K e L
              
                  F.[i] = Z. %*% P[,,i] %*% t(Z.) + H.
                  K.[,,i] = T. %*% P[,,i] %*% t(Z.) %*% solve(F.[i])
                  L.[,,i] = T. - matrix(K.[,,i]) %*% Z.
                  
              # Equacoes de atualizacao
                  
                  P[,,i+1] = T. %*% P[,,i] %*% t( T. - ( matrix(K.[,,i]) %*% Z. ) ) + R. %*% Q. %*% t(R.) 
                  a[,,i+1] = T. %*% matrix(a[,,i]) + matrix(K.[,,i]) %*% v[i]
                  
              # Estimador de y
                  
                  y.hat[i+1] = Z. %*% matrix(a[,,i+1])
              
            }

        
        # Coloca a serie do modelo estimado no grafico
            
            ts.plot(cbind(y,dropFirst(y.hat)),col=c("black","red"))
    
    # Suavizacao, recursao reversa, para estimar a variavel de estado com base em toda a informacao realizada
        
        # Notacoes auxiliares: alfa.hat, V, N, r
            
            alfa.hat = array(dim = c(2,1,t))
            V. = array(dim = c(2,2,t))
            
            r. = array(dim = c(2,1,t))
            N. = array(dim = c(2,2,t))
            
            r.[,,t] = matrix(c(0,0), nrow = 2, ncol = 1)
            N.[,,t] = matrix(rep(0,4), nrow = 2, ncol = 2)
            
        # Recursao para estimar y para i = t, ..., 1 por meio da suavizacao, conforme equacoes (4.44)
            
            # Declara o vetor y.hat.s
            
                y.hat.s = rep(0,t)
                
            for (i in t:1){
              
              if (i==1){ # Primeira observacao: Pela construcao da recursao reversa, fica fora do vetor
                
                # Auxiliares
                
                    r_0 = t(Z.) %*% solve(F.[i]) %*% v[i] + t(L.[,,i]) %*% matrix(r.[,,i])
                    N_0 = t(Z.) %*% solve(F.[i]) %*% Z. + t(L.[,,i]) %*% N.[,,i] %*% L.[,,i]
                
                # Estima o vetor de estado e sua variancia
                
                    alfa.hat[,,i] = a[,,i] + P[,,i] %*% matrix(r_0)
                    V.[,,i] = P[,,i] - P[,,i] %*% N_0 %*% P[,,i]
                
                # Estimador de y via suavizacao
                
                    y.hat.s[i] = Z. %*% matrix(alfa.hat[,,i])
                
              }
              else {
                
                # Auxiliares
                
                    r.[,,i-1] = t(Z.) %*% solve(F.[i]) %*% v[i] + t(L.[,,i]) %*% matrix(r.[,,i])
                    N.[,,i-1] = t(Z.) %*% solve(F.[i]) %*% Z. + t(L.[,,i]) %*% N.[,,i] %*% L.[,,i]
                    
                # Estima o vetor de estado e sua variancia
                
                    alfa.hat[,,i] = a[,,i] + P[,,i] %*% matrix(r.[,,i-1])
                    V.[,,i] = P[,,i] - P[,,i] %*% N.[,,i-1] %*% P[,,i]
                    
                # Estimador de y via suavizacao
                
                    y.hat.s[i] = Z. %*% matrix(alfa.hat[,,i])
                
              }
              
            }
        
        # Coloca a serie do modelo estimado (suavizado) no grafico
        ts.plot(cbind(y,(dropFirst(y.hat)),(y.hat.s)),col=c("black","red","blue"))
        ts.plot(cbind(y,(y.hat.s)),col=c("black","blue"))
        
# Plots mais ajeitados
        
        # Plot da serie original
        plot.ts(y, col = "darkgrey", xlab="",ylab = "GOOG",pch=3,cex=0.5,
                cex.lab=0.8,cex.axis=0.7)
        
        # Adiciona o modelo estimado
        lines(dropFirst(y.hat) , col = "red", lwd = 2, lty=2)
        
        # Legenda
        legend("topleft",leg = c("Preço de fechamento","Modelo estimado"),
               cex = 0.6,lty = c(1, 2), col = c("darkgrey","red"), pch=c(3,NA), bty = "y", horiz = T)
        
        # Plot do componente de nível
        plot.ts(alfa.hat[1,,], col = "darkgrey", xlab="",ylab = "")
        
        # Plot do componente de tendencia
        plot.ts(alfa.hat[2,,], col = "darkgrey", xlab="",ylab = "")
        
        # Plot do componente irregular (epsilon)
        plot.ts((y-y.hat.s), col = "darkgrey", xlab="",ylab = "")
        