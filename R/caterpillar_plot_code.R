


# chain plots (caterpillars)
par(mfrow = c(2,2))

pdf("my_output.pdf")
for(i in 1:ncol(chain_output[[1]])){
  my_mcmc <- sapply(
    chain_output,
    function(x){
      x[,i]
    }
  )
  
  my_range <- range(my_mcmc)
  
  plot(my_mcmc[,1], ylim = my_range, main = colnames(chain_output[[1]])[i],
       col = alpha('black', 0.5), type = 'l', ylab = 'value', xlab = "step")
  lines(my_mcmc[,2], col = alpha('red', 0.5))
  lines(my_mcmc[,3], col = alpha('green', 0.5))
  lines(my_mcmc[,4], col = alpha('purple', 0.5))
}
dev.off()

