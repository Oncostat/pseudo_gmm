data {
  int <lower=0> n; // nb of individuals
  int <lower=0> K; // nb of clusters
  real N; // nb of individuals
  int< lower = 0> np; // np of parameters
  matrix[n*K, np] X; // covariable matrix 
  vector[n*K] Y; // outcome variable
}

parameters {
vector[np] beta; // vector of parameters <lower = -5, upper = 5>
}


model {
  int loc1 = 0;
  int loc2 = 0;
  vector [K] yi; // outcome vector
  matrix [K, np] xi; // covariable matrix 
  vector [K] mui; // mean vector
  matrix[K,K] mui_dev; // derivation of the mean vector
  matrix[np,K] Di; // Di transpose  (!confusing notation!)
  vector[np] ui; // (1/n)*ui  (!confusing notation!)
  vector[np] U = rep_vector(0,np); // score vector
  matrix[np,np] C = rep_matrix(0, np, np); 
  matrix[np,np] Sigma; // empirical variance-covariance matrix


beta~ normal(0, sqrt(10)); //or normal(0, 1);
      
for(i in 1:n){
  
      // set start location for next xi
      loc1 = loc2+1;
      loc2 = loc1+K-1;
      
      
      yi = Y[loc1:loc2];
      xi = X[loc1:loc2,];
      
     mui = exp(-exp(xi*beta));
     mui_dev = diag_matrix(-exp(xi*beta))*diag_matrix(mui);

     Di = xi'*mui_dev;
     ui = (1/N)*Di*(yi-mui);
     C = C + ui*ui';
     U = U + ui;
        
        
}
Sigma = C - (1/N)*U*U';
        //pseudolikelihood for GMM :
        target += -0.5*U'/Sigma*U;
}

generated quantities {
  real loglik; 
  {
  int loc1 = 0;
  int loc2 = 0;
  vector [K] yi;
  matrix [K, np] xi; // covariable matrix variable
  vector [K] mui;
  matrix[K,K] mui_dev;
  matrix[np,K] Di;
  vector[np] ui;
  vector[np] U = rep_vector(0,np);
  matrix[np,np] C = rep_matrix(0, np, np);
  matrix[np,np] Sigma;

for(i in 1:n){

      // set start location for next xi
      loc1 = loc2+1;
      loc2 = loc1+K-1;
      
      
      yi = Y[loc1:loc2];
      xi = X[loc1:loc2,];
      
     mui = exp(-exp(xi*beta));
     mui_dev = diag_matrix(-exp(xi*beta))*diag_matrix(mui);

     Di = xi'*mui_dev;
     ui = (1/N)*Di*(yi-mui);
     C = C + ui*ui';
     U = U + ui;
}
    Sigma = C - (1/N)*U*U';
    loglik = -0.5*U'/Sigma*U;// added to the output

  }
}

