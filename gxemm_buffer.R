rm( list=ls() )
library(GxEMM)

savefile <- 'gxemm_buffer.Rdata'
if( file.exists( savefile ) ){
  load( savefile )
  for( kk in 1:2 ){

    print( nobs   <- sum( !is.na( (out['sig2g_1',kk,,]+out['sig2g_2',kk,,]) ) ) )

    print( mean( 
    (out['sig2g_1',kk,,]+out['sig2g_2',kk,,])/
    (out['sig2g_1',kk,,]+out['sig2g_2',kk,,]+out['sig2e_1',kk,,]+out['sig2e_2',kk,,])
    , na.rm=T ) )

    print( sd( 
    (out['sig2g_1',kk,,]+out['sig2g_2',kk,,])/
    (out['sig2g_1',kk,,]+out['sig2g_2',kk,,]+out['sig2e_1',kk,,]+out['sig2e_2',kk,,])
    , na.rm=T ) /sqrt( nobs ) )

  }
  stop()
}

 
sink( 'gxemm_buffer.Rout' )
N <- 1e3 # sample size
S <- 1e2 # number SNPs
threshes  <- seq(.5,.99,len=7)
maxit <- 1e2

sig2gs  <- c(.2,.5)

prev <- .16

out <- array( NA, dimnames=list( c( 'h2_1', 'h2_2', 'sig2g_hom','sig2g_1','sig2g_2','sig2e_1','sig2e_2'), 1:2, threshes, 1:maxit ), dim=c(7,length(sig2gs),length(threshes),maxit) )
for( j in 1:maxit )
  for( jj in 1:length(threshes) )
    for( kk in 1:2 )
{
  set.seed(1234+j+jj*1e3)  
  
  sig2g <- sig2gs[kk]
  G   <- scale( matrix( rbinom(N*S,2,.1), N, S ) )
  K   <- G %*% t(G) / S  
  betas	<- rnorm( S, sd=sqrt( sig2g/S ) )

  y0   <- as.numeric( G %*% betas + sqrt(1-sig2g) * rnorm(N)  )
  
  tau <- quantile(y0,1-prev)
  Z   <- as.numeric( y0 > tau )

  y   <- y0 + Z * (-3)

  out_free	<- GxEMM_HE( y, Z, K, cbind( Z, 1-Z ), gtype='free', etype='free' )
  print( c( out_free$h2, out_free$sig2g, out_free$sig2e ) )
  out[,kk,jj,j] <- c( out_free$h2, out_free$sig2g, out_free$sig2e )
  rm( out_free )

  print( apply( out, 1:2, mean, na.rm=T ) )
  save( out, file=savefile )

  if( sum( !is.na( (out['sig2g_1',kk,,]+out['sig2g_2',1,,]) ) )  == 20 &
      sum( !is.na( (out['sig2g_1',kk,,]+out['sig2g_2',2,,]) ) )  == 20   ) break

}
sink()
