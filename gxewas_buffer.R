rm( list=ls() )

N <- 1e4 # sample size
S <- 1e2 # number SNPs
h2<- .5
types <- c( 'add', 'thresh' )

pdf( 'gxewas_buffer.pdf', w=12, h=6 )
par( mfrow=c(2,4) )
par( mar=c(5,5,1,1) )
for( type in types ){
  set.seed(1234)  

  G   <- scale( matrix( rbinom(N*S,2,.1), N, S ) )
  K   <- G %*% t(G) / S  
  betas	<- rnorm( S, sd=sqrt( h2/S ) )
  y0   <- as.numeric( G %*% betas + sqrt(1-h2) * rnorm(N)  )

  tau <- quantile(y0,.8)
  Z   <- as.numeric( y0 > tau )

  if( type == 'add' ){
    y   <- y0 + Z * (-1)
  } else if( type == 'thresh' ){
    y <- y0
    y[y0>tau] <- tau
  }

  breaks <-  hist( y0, plot=F, breaks=51 )$breaks
  hist( y[Z==0], breaks=breaks, col=2, xlab='Post-Treatment Phenotypes', main='' )
  hist( y[Z==1], breaks=breaks, col=1, add=T ) 
  legend( 'topleft', fill=2:1, leg=c('Off Drug', 'On Drug'), bty='n' )

  b     <- apply( G, 2, function(snp) lm( y       ~ snp             )$coef[2] )
  b_gxe <- apply( G, 2, function(snp) lm( y       ~ snp + Z + snp*Z )$coef[4] )
  b0    <- apply( G, 2, function(snp) lm( y0      ~ snp             )$coef[2] )
  b_on  <- apply( G, 2, function(snp) lm( y[Z==1] ~ snp[Z==1]       )$coef[2] )
  b_off <- apply( G, 2, function(snp) lm( y[Z==0] ~ snp[Z==0]       )$coef[2] ) 

  lims <- max(abs(c(b,b_gxe,b0,b_on,b_off)))*c(-1,1)
  plot( b     , b_gxe , xlim=lims, ylim=lims, pch=16, xlab='Additive', ylab='GxDrug' )
  #plot( b0     , b_gxe , xlim=lims, ylim=lims, pch=16, xlab='Additive', ylab='GxDrug' )
  abline(h=0); abline(v=0)
  plot( b_off , b_on  , xlim=lims, ylim=lims, pch=16, xlab='Beta_off' , ylab='Beta_on' )
  abline(h=0); abline(v=0); abline(a=0,b=1)
  plot( b     , b0    , xlim=lims, ylim=lims, pch=16, xlab='Additive', ylab='Additive (Pre-Treat)' )
  abline(a=0,b=1)

  #plot( b0  , betas   , xlim=lims, ylim=lims, pch=16, xlab='Additive (Pre-Treat)', ylab='Truth' )
  #abline(a=0,b=1)

}
dev.off()
