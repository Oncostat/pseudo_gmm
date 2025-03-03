//Bayesian Generalized Method of Moments Model applied to (adjusted) RMST estimation
data {
  int <lower=0> n; // nb of individuals
  real N; // nb of individuals
  int< lower = 0> np; // np of parameters
  matrix[n, np] X; // covariable matrix 
  vector[n] Y; // outcome variable
  //vector[np] init; // init values
}

parameters {
vector[np] beta; // vector of parameters
}


model {
  real yi; // outcome value for one individual
  row_vector[np] xi; // covariable vector for one ind.
  real mui; // mean value
  //real mui_dev; // derivation of the mean vector
  vector[np] Di; // Di transpose  (!confusing notation!)
  vector[np] ui; // (1/n)*ui  (!confusing notation!)
  vector[np] U = rep_vector(0,np); // score vector
  matrix[np,np] C = rep_matrix(0, np, np); 
  matrix[np,np] Sigma; // empirical variance-covariance matrix

beta ~ normal(0, sqrt(10)); //normal(0, 10); //


      
for(i in 1:n){
  
      yi = Y[i];
      xi = X[i,];
      
     mui = xi*beta;
     //mui_dev = identity_matrix(K);

     //Di = xi'; //Di = xi'*mui_dev;
     ui = (1/N)*xi'*(yi-mui);
     C = C + ui*ui';
     U = U + ui;
}

Sigma = C - (1/N)*U*U';
        //pseudolikelihood for GMM :
        
target += -0.5*U'/Sigma*U;
//target += -0.5*mdivide_right_spd(U', Sigma)*U;
}

generated quantities {
  real loglik; 
  matrix[np,np] Sigma_n;
  matrix[np,np] C_n;
  {
  real yi; // outcome value for one individual
  row_vector[np] xi; // covariable vector for one ind.
  real mui; // mean value
  //real mui_dev; // derivation of the mean vector
  vector[np] Di; // Di transpose  (!confusing notation!)
  vector[np] ui; // (1/n)*ui  (!confusing notation!)
  vector[np] U = rep_vector(0,np); // score vector
  matrix[np,np] C = rep_matrix(0, np, np); 
  matrix[np,np] Sigma; // empirical variance-covariance matrix

for(i in 1:n){

      yi = Y[i];
      xi = X[i,];
      
      
     mui = xi*beta;
     //mui_dev = identity_matrix(K);

     //Di = xi'; //Di = xi'*mui_dev;
     ui = (1/N)*xi'*(yi-mui);
     C = C + ui*ui';
     U = U + ui;
}

    Sigma = C - (1/N)*U*U';
    C_n = C;
    Sigma_n = Sigma;
    //loglik =  -0.5*mdivide_right_spd(U', Sigma)*U; // added to the output
    loglik = -0.5*U'/Sigma*U;

  }
}

