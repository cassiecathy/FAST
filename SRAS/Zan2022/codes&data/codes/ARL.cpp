#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <queue>
using namespace Rcpp;
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//FUNCTION

// calculate the unit vector
mat SpatialSign (vec x) {
  double d = norm (x);
  if (d != 0) {
    return x / d;
  } else {
    return x;
  }
}

//
// generate data following different distributions
//

//
// generate data following certain distribution with fixed mean and covariance
mat Generate (string Distribution, int n, vec mu, mat sigma){
  // 
  // Distribution: can be chosen from {norm, t, gamma}
  // n:# of data;
  // mu: mean vector;
  // sigma: covariance matrix
  //
  if (Distribution == "norm") {
    return rmvnorm(n, mu, sigma);
  } else if (Distribution == "t") {
    double df = 3.0;
    sigma = sigma * (df - 2) / df;
    return rmvt(n, mu, sigma, df);
  }
  else {
    cout<<"WRONG DISTRIBUTION";
    return 0;
  }
}

// zca decomposition
mat zca(mat sigma) {
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, sigma);
  int n = eigval.size();
  mat diagmat(n, n);
  diagmat.fill(0);
  for (int i = 0; i < n; i++){
    diagmat(i,i) = 1 / sqrt(eigval(i));
  }
  mat M_zca = eigvec * diagmat * eigvec.t();
  return M_zca;
}


// Determine which r variables are chosen observable at time t, update of layout
IntegerVector LAYOUT (vec v, int n) {
  typedef pair<double, int> Elt;
  priority_queue< Elt, vector<Elt>, greater<Elt> > pq;
  vector<int> result;
  
  for (int i = 0; (unsigned) i != v.size(); ++i) {
    if (pq.size() < (unsigned) n) {
      pq.push (Elt(abs(v[i]), i));
    }
    else {
      Elt elt = Elt(abs(v[i]), i);
      if (pq.top() < elt) {
        pq.pop();
        pq.push(elt);
      }
    }
  }
  
  result.reserve (pq.size());
  while (!pq.empty()) {
    result.push_back (pq.top().second);
    pq.pop();
  }
  
  return wrap(result);
}

// NW estimation function
vec mNWcpp (NumericVector x, NumericMatrix X, NumericMatrix Y, NumericVector h) {
  //
  // x: a single input observation, x.size = # of var = X.ncol
  // X: historical observations, X.nrow = # of observations, X.ncol = # of input var
  // Y: historical observations, Y.nrow = X.nrow = # of observations, Y.ncol = # of output responses
  // h : bandwidth for each Xi, h.size = # of var = X.ncol = x.size
  //
  if (X.nrow() != Y.nrow() || X.ncol() != x.size() || X.ncol() != h.size()){
    cout<<"something goes wrong"<<endl;
  }
  int number_loop = X.nrow();
  int number_var = X.ncol();
  
  NumericVector Kx(number_loop);
  Kx.fill(1);
  
  for (int i = 0; i < number_loop; ++i){
    for (int j = 0; j < number_var; ++j){
      // Kx(i) = Kx(i) * R::dnorm((x(j)-X(i,j)) / h(j), 0.0, 1.0, false);
      Kx(i) = Kx(i) * exp(- pow((x(j) - X(i,j)) / h(j), 2) / 2);
    }
  }
  
  double row_sum = sum(Kx);
  if (row_sum == 0){
    row_sum = 1e-300;
  }
  NumericVector W = Kx / row_sum;
  
  
  vec weights = as<vec>(W);
  mat Ymat = as<mat>(Y);
  
  vec y = Ymat.t() * weights;
  
  return y;
  // y.size = # of output responses
}

vec Lsquare(vec Rf){
  int n = Rf.size();
  vec L(n);
  for (int i = 0; i < n; i++){
    L(i) = pow(Rf(i), 2); 
  }
  return L;
}

// monitoring process: function for calculating average run length (both IC and OC)


// [[Rcpp::export]]
mat ARL(const int q, const int m0, const int n, const double h,
        const int p, const int r, const double lambda, const double L,
        string Distribution, vec mu, mat sigma, double shift, double delta){
  //
  // q:# of experiment times
  // m0: # of offline observations
  // n:# of periods
  // p:# of monitoring variables
  // r:# of observable variables
  // L: control limit
  // Distribution: can be chosen from {norm, t, gamma}
  // mu: in-control mean vector
  // sigma: in-control covariance matrix
  // shift: out-of-control mean shift
  // delta: complementary penalty
  //
  
  vec Arl(q);
  Arl.fill(n);
  
  // bandwidth vector
  NumericVector H(r);
  H.fill(h);
  
  // // create a matrix to record the observed frequency of each variable
  // mat Sum_Freq(p, n);
  // Sum_Freq.fill(0);
  
  // create a matrix to record Rfs and Ws
  mat Rf_MX(p, n);
  mat Rf_MX_off(p, m0);
  mat L_Rf(p, n);
  mat W(p, n);
  mat Z(p, n);
  
  // create a vec to record charting stat Q
  vec Q(n);
  
  for (int k = 0; k < q; k++){
    // generate m0 historical IC observations randomly according to certain distribution
    // calculate initial estimators
    mat X_offline = Generate(Distribution, m0, mu, sigma);
    mat mu_0 = mean(X_offline);
    mat cov_0 = cov(X_offline);
    mat M_0 = zca(cov_0);
    
    // "online monitoring procedure"
    
    // statistics setting
    double AVR_pow_Rf = 0;
    mat U_MX_off(p, m0);
    for (int j = 0; j < m0; j++){
      U_MX_off.fill(0);
      for (int i = 0; i < m0; i++){
        U_MX_off.col(i) = SpatialSign(M_0 * (X_offline.row(j) - X_offline.row(i)).t());
      }
      Rf_MX_off.col(j) = mean(U_MX_off, 1);
      AVR_pow_Rf += pow(norm(Rf_MX_off.col(j)), 2);
    }
    AVR_pow_Rf = AVR_pow_Rf / m0;
    Rf_MX.fill(0);
    L_Rf.fill(0);
    W.fill(0);
    Z.fill(0);
    Q.fill(0);
    
    // set the OC mu_1 to make the one with shift chosen randomly each time
    vec mu_1 = mu;
    IntegerVector pool = seq(0, p - 1);
    random_shuffle(pool.begin(), pool.end());
    int v = pool[0];
    mu_1[v] += shift;
    
    // generate data randomly, following multivariate normal distribution with shift
    mat X_online = Generate(Distribution, n, mu_1, sigma);
    mat X_observe = X_online;
    
    //generate initial layout randomly
    random_shuffle(pool.begin(), pool.end());
    IntegerVector layout = pool[Range(0, r - 1)]; // layout indicates which r varibles are observed
    
    
    // online data statistics monitoring process
    for (int t = 0; t < n; t++){
      
      mat X_off(m0,r);
      mat Y_off(m0,p-r);
      NumericVector x_ob(r);
      
      int a = 0;
      int b = 0;
      int c = 0;
      
      for (int i = 0; i < p; i++){
        bool unobserved = find(layout.begin(), layout.end(), i) == layout.end();
        
        // update unobserved X replaced with estimated mean
        if (unobserved){
          Y_off.col(b) = X_offline.col(i);
          ++b;
        }
        else{
          X_off.col(a) = X_offline.col(i);
          x_ob(a) = X_online(t,i);
          ++a;
        }
      }
      
      vec y = mNWcpp (x_ob, wrap(X_off), wrap(Y_off), H);
      
      for (int i = 0; i < p; i++){
        
        // update unobserved X replaced with estimated mean
        if (find(layout.begin(), layout.end(), i) == layout.end()){
          if (NumericVector::is_na(y(c))){
            cout<<"something is na during NW estimation";
          }
          X_observe(t,i) = y(c);
          c++;
        }
      }
      
      
      // calculate Rf_MX
      mat U_MX_online(p, m0);
      for (int i = 0; i < m0; i++){
        U_MX_online.col(i) = SpatialSign(M_0 * (X_observe.row(t) - X_offline.row(i)).t());
      }
      
      Rf_MX.col(t) = mean(U_MX_online, 1);
      
      // L change
      L_Rf.col(t) = Lsquare(Rf_MX.col(t));
      
      // penalty to L
      for (int i = 0; i < p; i++){
        if (find(layout.begin(), layout.end(), i) == layout.end()){
          L_Rf(i,t) += delta;
        }
      }
      
      // update AVR_pow_Rf and ewma sequence statistic
      if (t == 0){
        Z.col(t) = lambda * L_Rf.col(t);
        W.col(t) = lambda * Rf_MX.col(t);
      }else{
        Z.col(t) = lambda * L_Rf.col(t) + (1 - lambda) * Z.col(t - 1);
        W.col(t) = lambda * Rf_MX.col(t) + (1 - lambda) * W.col(t - 1);
      }
      
      // update layout, i.e. find the most possible out-of-control variables, say, the top r largest element of w
      layout = LAYOUT(Z.col(t), r);
      
      // update Q at time t
      Q(t) = (2 - lambda) * p / (lambda * AVR_pow_Rf) * pow(norm(W.col(t)), 2);
      
      // judge whether in control
      if ((Q(t) > L) || t == n - 1){
        Arl(k) = t + 1;
        
        // cout << X_observe.rows(tau,t);
        
        // set the rest of freq in the following periods be the last freq
        // if (t < n - 1){
        //   for (int i = 0; i < p; i++){
        //     double converge = Freq(i, t);
        //     for (int j = t + 1; j < n; j++){
        //       Freq(i, j) = converge;
        //       Sum_Freq(i, j) += Freq(i, j);
        //     }
        //   }
        // }
        break;
      }
    }//run n periods
    
    // cout<<mu_0;
    // cout<<M_0;
    
  }//repeat q times
  
  
  vec RL(2);
  RL(0) = mean(Arl);
  RL(1) = stddev(Arl) / sqrt(q);
  // cout<<RL;
  return RL;
}