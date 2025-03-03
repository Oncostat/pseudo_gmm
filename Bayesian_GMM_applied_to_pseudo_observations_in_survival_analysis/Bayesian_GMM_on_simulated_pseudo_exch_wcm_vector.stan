data {
  int <lower=0> n; // nb of individuals
  int <lower=0> K; // nb of clusters
  real N; // nb of individuals
  int< lower = 0> np; // np of parameters
  matrix[n*K, np] X; // covariable matrix variable
  vector[n*K] Y; // outcome variable
}

parameters {
vector[np] beta; // vector of parameters
}


model {
  int loc1 = 0;
  int loc2 = 0;
  vector [K] yi;
  matrix [K, np] xi; // covariable matrix variable
  vector [K] ui;
  matrix[K,K] fui_dev;
  matrix[K,K] M2;
  matrix[np,K] Di;
  vector[2*np] gi;
  vector[2*np] U = rep_vector(0,2*np);
  matrix[2*np,2*np] C = rep_matrix(0, 2*np, 2*np);
  matrix[2*np,2*np] Sigma;


//priors:
beta~ normal(0, sqrt(10)); //normal(0, 1);


M2[1,1] = 0.0; M2[1,2] = 1.0; M2[1,3] = 1.0; M2[1,4] = 1.0; M2[1,5] = 1.0;
M2[2,1] = 1.0; M2[2,2] = 0.0; M2[2,3] = 1.0; M2[2,4] = 1.0; M2[2,5] = 1.0;
M2[3,1] = 1.0; M2[3,2] = 1.0; M2[3,3] = 0.0; M2[3,4] = 1.0; M2[3,5] = 1.0;
M2[4,1] = 1.0; M2[4,2] = 1.0; M2[4,3] = 1.0; M2[4,4] = 0.0; M2[4,5] = 1.0;
M2[5,1] = 1.0; M2[5,2] = 1.0; M2[5,3] = 1.0; M2[5,4] = 1.0; M2[5,5] = 0.0;

for(i in 1:n){
  
      // set start location for next xi
      loc1 = loc2+1;
      loc2 = loc1+K-1;
      
      
      yi = Y[loc1:loc2];
      xi = X[loc1:loc2,];
      
     ui = exp(-exp(xi*beta));
     fui_dev = diag_matrix(-exp(xi*beta))*diag_matrix(ui);

     Di = xi'*fui_dev;
     gi[1:np] = (1/N)*Di*(yi-ui);
     gi[(np+1):2*np] = (1/N)*Di*M2*(yi-ui);
     C = C + gi*gi';
     U = U + gi;
        
        
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
  vector [K] ui;
  matrix[K,K] fui_dev;
  matrix[np,K] Di;
  matrix[K,K] M2;
  vector[2*np] gi;
  vector[2*np] U = rep_vector(0,2*np);
  matrix[2*np,2*np] C = rep_matrix(0, 2*np, 2*np);
  matrix[2*np,2*np] Sigma;

M2[1,1] = 0.0; M2[1,2] = 1.0; M2[1,3] = 1.0; M2[1,4] = 1.0; M2[1,5] = 1.0;
M2[2,1] = 1.0; M2[2,2] = 0.0; M2[2,3] = 1.0; M2[2,4] = 1.0; M2[2,5] = 1.0;
M2[3,1] = 1.0; M2[3,2] = 1.0; M2[3,3] = 0.0; M2[3,4] = 1.0; M2[3,5] = 1.0;
M2[4,1] = 1.0; M2[4,2] = 1.0; M2[4,3] = 1.0; M2[4,4] = 0.0; M2[4,5] = 1.0;
M2[5,1] = 1.0; M2[5,2] = 1.0; M2[5,3] = 1.0; M2[5,4] = 1.0; M2[5,5] = 0.0;

for(i in 1:n){

      // set start location for next xi
      loc1 = loc2+1;
      loc2 = loc1+K-1;
      
      
      yi = Y[loc1:loc2];
      xi = X[loc1:loc2,];
      
     ui = exp(-exp(xi*beta));
     fui_dev = diag_matrix(-exp(xi*beta))*diag_matrix(ui);

     Di = xi'*fui_dev;
     gi[1:np] = (1/N)*Di*(yi-ui);
     gi[(np+1):2*np] = (1/N)*Di*M2*(yi-ui);
     C = C + gi*gi';
     U = U + gi;
}

    Sigma = C - (1/N)*U*U';
    loglik = -0.5*U'/Sigma*U; // added to the output
  }
}
