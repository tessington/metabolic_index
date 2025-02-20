#include <TMB.hpp>
  
  // Generate covariance matrix
  // If positive, then factor decomposition
  // If zero, then diagonal matrix
  // If negative, then both factor decomposition and diagonal components
  template<class Type>
  matrix<Type> cov_matrix( vector<Type> L_val, Type min_var, int n_rows){
    // Define temporary objects
    matrix<Type> L_rc(n_rows, n_rows);
    matrix<Type> Cov_rr(n_rows, n_rows);
    matrix<Type> Return_rr(n_rows, n_rows);
    Cov_rr.setZero();
    L_rc.setZero();
    // Loadings matrix with zero upper-diagonal
    int Count = 0;
    
      for(int r=0; r<n_rows; r++){
        for(int c=0; c<n_rows; c++){
          if(r>=c){
            L_rc(r,c) = L_val(Count);
            Count++;
          }else{
            L_rc(r,c) = 0.0;
          }
        }}
    
    
    // Additive constant on diagonal
    for(int r=0; r<n_rows; r++){
      Cov_rr(r,r) += min_var;  // Neccesary to prevent crashes during innner optimizer when using SR data
    }
    // Combine and return
    Cov_rr += L_rc * L_rc.transpose();
     Return_rr = Cov_rr;
    return Return_rr;
  }
  
  // main TMB 
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
   
    
    // Data 
    
    DATA_IMATRIX( PC_gz );
    DATA_IVECTOR( g_i );
    DATA_VECTOR( invtemp );
    DATA_VECTOR( logW );
    DATA_IVECTOR( taxa_id );
    DATA_VECTOR( minuslogpo2 );
    DATA_IVECTOR( spc_in_PCgz );
  
    
    // Parameters 
    PARAMETER_VECTOR( alpha_j );
    PARAMETER_VECTOR( L_z );
    PARAMETER_VECTOR( log_lambda ); // log-multiplier for process-error covariance for different taxonomic levels
    PARAMETER_MATRIX( beta_gj );
    PARAMETER( logsigma );

    
    
    // Derived data and parameters
    Type sigma = exp(logsigma);
    int n_j = beta_gj.row(0).size();
    int n_i = spc_in_PCgz.size();
    int n_g = PC_gz.col(0).size();
    int n_d = minuslogpo2.size();
    int n_l = log_lambda.size();
    
  
    vector<Type> Parent_j( n_j );
    vector<Type> Prediction_j( n_j );
    vector<Type> Deviation_j( n_j );
    matrix<Type> tmpCov_jj( n_j, n_j );
    vector<Type> V(n_i);
    vector<Type> n_pow(n_i);
    vector<Type> Eo(n_i);
    matrix<Type> spc_ij(n_i, n_j);
    vector<Type> lambda(n_l);
    
    
    // to make life easier, extract the species-level parameters into a matrix, and convert logAo to Ao
    
    for (int i = 0; i <n_i; i++) {
      for (int j = 0; j<n_j; j++) {
        if (j == 0) {
          spc_ij(i,j) = exp(beta_gj(spc_in_PCgz( i ), j ));
        }
        if (j >0) {
          spc_ij(i,j) = beta_gj(spc_in_PCgz( i ), j );
        }
      }
    }
    // make vector of lambdas
    for(int l = 0; l < n_l; l++) {
       lambda( l ) =  exp(log_lambda( l ));
    }
    
    V = spc_ij.col(0);
    n_pow = spc_ij.col(1);
    Eo = spc_ij.col(2);
    vector<Type> mu( n_d );
    
    // Objective function
    vector<Type> jnll_comp( 2 );
    jnll_comp.setZero();
      
    using namespace density;
    
    // Process covariance
    Type min_var = 0.001;
    matrix<Type> Cov_jj( n_j, n_j );
    Cov_jj = cov_matrix(L_z, min_var, n_j);
    
    // Probability of random effects above species level
   Type covmult;
   
    for( int g=0; g<n_g; g++ ){
      int Child_num = PC_gz(g,1);
      int Parent_row = PC_gz(g,0);
      int lambda_num = Child_num - 1;
      for( int j=0; j<n_j; j++ ){
        if( PC_gz(g,1)==0 ) Parent_j(j) = alpha_j(j);
        if( PC_gz(g,1)>=1 ) Parent_j(j) = beta_gj(Parent_row,j);
        Prediction_j(j) = Parent_j(j);
      }
      for( int j=0; j<n_j; j++ ){
        Deviation_j(j) = beta_gj(g,j) - Prediction_j(j);
      }
      if( Child_num==0 ) covmult = 1;
      if( Child_num>=1 ) covmult = lambda( lambda_num );
      tmpCov_jj = Cov_jj * covmult;
        jnll_comp(0) += MVNORM( tmpCov_jj )( Deviation_j );
      }
    
    
    // Probability of the data
    for( int id=0; id<n_d; id++){
     mu( id ) =  Eo( taxa_id( id ) ) * invtemp( id ) + n_pow( taxa_id( id ) ) * logW( id  ) - log(V( taxa_id( id ) )) ;
    }
    
    jnll_comp( 1 ) = -sum( dnorm( minuslogpo2, mu, sigma, true) );
    Type jnll = jnll_comp.sum();
    
    // Return jnll
    // Add Reports
    ADREPORT(spc_ij);
    return jnll;
  }
  
