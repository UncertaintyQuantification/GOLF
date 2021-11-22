
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 
 
// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <cmath> 
#include "ctools.h"
#include <random>
#include <cstdlib>


// [[Rcpp::depends(RcppEigen)]] 


//#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
// <Spectra/MatOp/DenseSymMatProd.h> is implicitly included
//#include <iostream>

using namespace Rcpp;
using namespace std;
using namespace Eigen; 
 
using namespace Spectra;
 

//July 16, 2018
////Construct_W0_matern_5_2 
// [[Rcpp::export]] 
MatrixXd Construct_W0_matern_5_2(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(3,3); 
  //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  W0(0,0)=sigma2; 
  W0(0,2)=W0(2,0)=-sigma2*pow(lambda,2.0)/3.0; 
  W0(1,1)=sigma2*pow(lambda,2.0)/3.0; 
  W0(2,2)=sigma2*pow(lambda,4.0); 

  return W0; 
} 

//July 16, 2018
////Construct_W0_exp 
// [[Rcpp::export]] 
MatrixXd Construct_W0_exp(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(1,1); 

  W0(0,0)=sigma2; 
  
  return W0; 
} 


 
 ////Construct_G_matern_5_2 
 // [[Rcpp::export]] 
 List Construct_G_matern_5_2(const Eigen::VectorXd &delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
   int num_obs=delta_x.size()+1; 
   //int num_dim=lambda.size();  
   List GG(num_obs);  
   GG[0]=Eigen::MatrixXd::Zero(3,3); 
   
   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
   
   // num_dim list, each is 3(num_obs)\times 3 list 
  // for(int i_GG=0;i_GG<num_dim;i_GG++){ 
  //   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);  //the first row has all zeros  
     for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
       int j_GG_1=j_GG+1;    
       d(0,0)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)+2*lambda*delta_x[j_GG]+2; 
       d(1,0)=-pow(lambda,3.0)*pow(delta_x[j_GG],2.0); 
       d(2,0)=pow(lambda,4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda,3.0)*delta_x[j_GG]; 
       d(0,1)=2*(lambda*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
       d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-lambda*delta_x[j_GG]-1); 
       d(2,1)=2*(pow(lambda,3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda,2.0)*delta_x[j_GG]); 
       d(0,2)=pow(delta_x[j_GG],2); 
       d(1,2)=2*delta_x[j_GG]-lambda*pow(delta_x[j_GG],2.0); 
       d(2,2)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-4*lambda*delta_x[j_GG]+2;     
       d=exp(-lambda*delta_x[j_GG])/2.0*d;
       GG[j_GG_1]=d; 
     } 
    // GG[i_GG]=d; 
   //} 
   return GG; 
} 

////Construct_G_exp
// [[Rcpp::export]] 
List Construct_G_exp(const Eigen::VectorXd &delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=lambda.size();  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(1,1); 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    d(0,0)=exp(-delta_x[j_GG]*lambda);
    GG[j_GG+1]=d; 
  }

  return GG;
}


////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_matern_5_2(double sigma2,const Eigen::VectorXd& delta_x, double lambda, const MatrixXd& W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  
 // List Wi(num_obs); 
  //for(int i_Wi=0;i_Wi<num_dim;i_Wi++){ 
    //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);   
    double  lambda_delta_x;
    double exp_neg_2_lambda_delta_x;
    int  j_Wi_1;
    for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
      j_Wi_1= j_Wi+1; 
      lambda_delta_x=lambda*delta_x[j_Wi];  //close and jump then it is... 
      exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
      
      d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
      d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],4.0)/2.0; 
      d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
      d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
      d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
      d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
      d=d*(4*sigma2*pow(lambda,5.0)/3.0); 
      Wi[j_Wi_1]=d; 
      
    //} 
  } 
  return Wi; 
}

////Construct_W_exp
// [[Rcpp::export]] 
List Construct_W_exp(double sigma2, const Eigen::VectorXd& delta_x, double lambda, const MatrixXd& W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    d(0,0)=1-exp(-2*delta_x[j_Wi]*lambda);
    Wi[j_Wi+1]=d;
  }
  
  return Wi;
}


////Get_Q_K  
// [[Rcpp::export]] 
List Get_Q_K(const List& GG,const List&  W,const Eigen::MatrixXd& C0,const double VV){ 

   int n=GG.size();
   int k=C0.rows();
   
   Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
   Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
   Eigen::MatrixXd C=C0;
    
   Eigen::MatrixXd GG_matrix;
   Eigen::MatrixXd W_matrix;
   
   Eigen::MatrixXd RR;
   
      
   // num_dim list, each is 3(num_obs)\times 3 list 
   for(int t=0;t<n;t++){ 
     GG_matrix=GG[t];
     W_matrix=W[t];
     RR=GG_matrix*C*GG_matrix.transpose()+W_matrix;
     //Q[t]=RR(0,0);
     Q[t]=RR(0,0)+VV;
     K.row(t)=RR.col(0).transpose()/Q[t];
     C=RR-RR.col(0)*RR.row(0)/Q[t];
   }

   List return_list;
   return_list.push_back(Q);
   return_list.push_back(K);
   
   return return_list;
}







////Get_C_R_K_pred, C and R is for smoothing, K is for filtering the mean  
// let me also output Q here for determinant
// [[Rcpp::export]] 
List Get_C_R_K_Q_pred(const VectorXi& index, const List& GG,const List&  W,const Eigen::MatrixXd& C0, double VV){ 
  
  
  //index is a sequence where 0 means NA and 1 means observations
  int n=GG.size();//total number
  int k=C0.rows();
  
  //Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
  List C(n+1);
  C[0]=C0;
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd W_matrix;
  Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  
  Eigen::MatrixXd RR;
  
  List R(n);
  
  Eigen::MatrixXd C_cur=C[0];
  //int index_num=0;
  
  for(int t=0;t<n;t++){
    
    if(index[t]==1){
      //index_num=index_num+1;
      GG_matrix=GG[t];
      W_matrix=W[t];
      RR=GG_matrix*C_cur*GG_matrix.transpose()+W_matrix;
      R[t]=RR;
      
      //Q[t]=RR[t](0,0);
      Q[t]=RR(0,0)+VV;
      
      
      K.row(t)=RR.col(0).transpose()/Q[t];
      C[t+1]=RR-RR.col(0)*RR.row(0)/Q[t];
      C_cur=C[t+1];
      
    }else{
      GG_matrix=GG[t];
      W_matrix=W[t];
      
      R[t]=GG_matrix*C_cur*GG_matrix.transpose()+W_matrix;
      C[t+1]=C_cur=R[t];   
      
    }
  }
  
  List return_list;
  return_list.push_back(C);
  return_list.push_back(R);
  return_list.push_back(K);
  return_list.push_back(Q);
  
  return return_list;
}


////Get_Y_minus_a_1_scaled_matrix_2d  
// [[Rcpp::export]] 
List Get_m_a_pred(const VectorXi& index, const Eigen::VectorXd& output_vec,const List& GG,const Eigen::MatrixXd& K){
  
  //output_KF
  int n=GG.size();//total number
  
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  int k=K.cols();
  
  List m(n);  //m should have n+1 item but the first one is (0,0,0) so it's omitted.
  List a(n);
  Eigen::VectorXd a_vec;
  Eigen::VectorXd m_cur=Eigen::VectorXd::Zero(k); //3\times 1
  
  // =Eigen::MatrixXd::Zero(k,n2);
  
  //Eigen::MatrixXd Y_minus_a_1_scaled_matrix=Eigen::MatrixXd::Zero(n1,n2); 
  //Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::MatrixXd GG_matrix;
  
  //int index_num=0;
  int index_count=0;
  for(int t=0;t<n;t++){
    if(index[t]==1){
      //index_num=index_num+1;
      GG_matrix=GG[t];
      a_vec=GG_matrix*m_cur;
      a[t]=a_vec;
      m[t]=a_vec+K.row(t).transpose()*(output_vec[index_count]-a_vec[0]);
      m_cur=m[t];
      index_count=index_count+1;
    }else{
      GG_matrix=GG[t];
      m[t]=a[t]=GG_matrix*m_cur;
      m_cur=m[t];
    }
  }
  
  List return_list;
  
  return_list.push_back(m);
  return_list.push_back(a);
  
  return return_list;
}


////Kalman_smoother, the index is 0 and 1 where 0 means missing value
// [[Rcpp::export]] 
List Get_S_KK(const VectorXi& index,const List& GG,const List& C,const List& R){
  
  int n=GG.size();
  //int k=3; //matern_5_2
  
  List S(n);
  List KK(n-1);
  
  S[n-1]= C[n];
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd C_matrix;
  Eigen::MatrixXd R_matrix;
  Eigen::MatrixXd KK_matrix;
  
  Eigen::MatrixXd S_cur=C[n];
  for(int t=n-2;t>=0;t--){
    GG_matrix=GG[t+1];
    R_matrix=R[t+1];
    C_matrix=C[t+1];  // C should have n+1 items
    
    LLT<MatrixXd> lltOfR(R_matrix);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
    
    KK_matrix= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(GG_matrix*C_matrix))); 
    KK[t]=KK_matrix;
    
    S[t]=C_matrix-KK_matrix.transpose()*(R_matrix-S_cur)*KK_matrix;
    S_cur=S[t];
  }
  
  List return_list;
  
  return_list.push_back(S);
  return_list.push_back(KK);
  
  return return_list;
  
}

// [[Rcpp::export]] 
MatrixXd Get_s_1st(const List& m,const List& a,const List& C,const List& KK){
  
  int n=C.size()-1;
  //int k=3; //matern_5_2
  
  Eigen::VectorXd s=m[n-1];
  //int n_2=s.cols();
  VectorXd s_1st=Eigen::VectorXd::Zero(n);
  MatrixXd KK_matrix;
  s_1st[n-1]=s[0];
  VectorXd a_vec;
  VectorXd m_vec;
  for(int t=n-2;t>=0;t--){
    KK_matrix=KK[t];
    a_vec=a[t+1];
    m_vec=m[t];
    s=m_vec+KK_matrix.transpose()*(s-a_vec);
    s_1st[t]=s[0];
  }
  
  return s_1st; 
}










 






   
// Nov 6, 2019 this function output the L^{-1}y given the constructed matrix and  Q K

////Get_L_inv_y
// [[Rcpp::export]] 
VectorXd Get_L_inv_y(const List& GG,const double VV,const Eigen::VectorXd& Q, const Eigen::MatrixXd& K,  const Eigen::VectorXd& output){

  //Eigen::VectorXd Q=Q_K[0];
  //Eigen::MatrixXd K=Q_K[1];
  
  //double log_det_R=Q.array().log().sum();
  //List return_vec;
  //return_vec.push_back(log_det_R);
  
  //return_list.push_back(log_det_R);
    int n=GG.size();
  //int k=C0.rows();
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(n);
  
  Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    Y_minus_a_1[t]=(output[t]-a[0]);
    
    // Y_minus_a_1_scaled_vec[t]=(output[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(output[t]-a[0]);
  }
  
  ///sigma2_S2=(FRFt+sigma_2*eta)  ####using this will lead to |R|
  
  Eigen::VectorXd res=(Y_minus_a_1.array()/Q.array().sqrt()).matrix();
  return res;
  //double S2=(Y_minus_a_1.array()*Y_minus_a_1.array()/Q.array()).sum(); 

}




// [[Rcpp::export]] 
List Kalman_smoother_mean_sample(const Eigen::MatrixXd& W0,const List& GG,const List& W,const double VV,
                                 const VectorXi& index_obs,const VectorXd& output){
  
  // //int n1=output_KF.rows();
  // //int n2=output_KF.cols();
  // 
  // double gamma=1.0/exp(param[0]);
  // //double lambda=sqrt(5.0)/gamma;
  // 
  // 
  // double VV=0;
  // if(have_noise){
  //   VV=exp(param[1]); 
  // }
  int n=output.size();

  
  
  List C_R_K_Q=Get_C_R_K_Q_pred(index_obs, GG,W, W0, VV);
  
  
  List m_a=Get_m_a_pred(index_obs, output, GG,C_R_K_Q[2]);
  
  List S_KK=Get_S_KK(index_obs, GG,C_R_K_Q[0],C_R_K_Q[1]);
  
  //posterior mean
  MatrixXd s_1st=Get_s_1st(m_a[0],m_a[1], C_R_K_Q[0],S_KK[1]);
  
  //new
  //posterior sample
  VectorXd sample_theta_1_record=Eigen::VectorXd::Zero(n);;
  Eigen::VectorXd theta_sample_here;
  
  //Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd theta_minus_a; 
  
  Eigen::VectorXd random_norm=Eigen::VectorXd::Zero(3);
  for(int t=0;t<3;t++){
    random_norm(t)= R::rnorm(0,1.0);
  }
  
//the first sample
  
  List C=C_R_K_Q[0];
  List m=m_a[0];

  Eigen::VectorXd theta_sample;
  
  Eigen::MatrixXd C_here=C[n]; //it looks C has n+1 terms but m only have m terms?
  Eigen::VectorXd m_here=m[n-1];
  
  LLT<MatrixXd> lltOfC_here(C_here);    // compute the cholesky decomposition of R called lltofR
  MatrixXd L_C_here = lltOfC_here.matrixL();   //retrieve factor L  in the decomposition
  
  theta_sample_here=m_here+L_C_here*random_norm;
  sample_theta_1_record[n-1]=theta_sample_here[0];
  
  //done for sample 1
  List S_list=S_KK[0];
  List KK_list=S_KK[1];
  
  List a=m_a[1];
  VectorXd a_here;
  //VectorXd m_here;
  
  VectorXd h;
  MatrixXd KK_matrix;
  MatrixXd S_matrix;
  
    /*
    for(int t=n-2;t>=0;t--){
      KK_matrix=KK[t];
      a_vec=a[t+1];
      m_vec=m[t];
      s=m_vec+KK_matrix.transpose()*(s-a_vec);
      s_1st[t]=s[0];
    }
    */
  for(int t=n-2;t>=0;t--){
    
    for(int t_norm=0;t_norm<3;t_norm++){
      random_norm(t_norm)= R::rnorm(0,1.0);
    }
    
    KK_matrix=KK_list[t];
    a_here=a[t+1];
    m_here=m[t];
    h=m_here+KK_matrix.transpose()*(theta_sample_here-a_here);
    
    S_matrix=S_list[t];
    
    LLT<MatrixXd> lltOfS(S_matrix);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L_S = lltOfS.matrixL();   //retrieve factor L  in the decomposition
    
    theta_sample_here=h+L_S*random_norm;
    sample_theta_1_record[t]=theta_sample_here[0];
  }
  
  List return_list;
  /*
   List S= S_KK[0];
   MatrixXd S_matrix;
   VectorXd Var=VectorXd::Zero(n);
   for(int t=0;t<n;t++){
   S_matrix=S[t];
   Var[t]=S_matrix(0,0);
   }
   //Var=((Var.array()+VV)*sigma_2_hat).matrix();
   Var=((Var.array()+VV)).matrix();
   */
  
  return_list.push_back(s_1st);
  return_list.push_back(sample_theta_1_record);
  
  //let me also get Q out 
  return_list.push_back(C_R_K_Q[3]);
  
  
  /*
   List return_list;
   
   return_list.push_back(W);
   */
  return return_list;
  
}


// [[Rcpp::export]]
bool Accept_proposal(double r){
  if(r>=1){
    return true;
  }else{
    
    //construct a trivial random generator engine from a time-based seed:
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    //std::default_random_engine generator (seed);
    
    //std::uniform_real_distribution<> dis(0.0, 1.0);
    
    //double u=dis(generator);
    
    double u= R::runif(0.0,1.0);
    if(u<r){
      return true;
    }else{
      return false;
    }
  }
  
}



//  Nov 25 separate the sample and mean
// [[Rcpp::export]] 
List Kalman_smoother_mean_sample_known_C_R_K_Q_matern_5_2(const List& C_R_K_Q, const Eigen::MatrixXd& W0,const List& GG,const List& W,const double VV,
                                               const VectorXi& index_obs,const VectorXd& output,  const bool need_sample, const bool need_mean,const bool need_var){
  
  // //int n1=output_KF.rows();
  // //int n2=output_KF.cols();
  // 
  // double gamma=1.0/exp(param[0]);
  // //double lambda=sqrt(5.0)/gamma;
  // 
  // 
  // double VV=0;
  // if(have_noise){
  //   VV=exp(param[1]); 
  // }
  int n=output.size();
  
  
  //List C_R_K_Q=Get_C_R_K_Q_pred(index_obs, GG,W, W0, VV);
  
  
  List m_a=Get_m_a_pred(index_obs, output, GG,C_R_K_Q[2]);
  
  List S_KK=Get_S_KK(index_obs, GG,C_R_K_Q[0],C_R_K_Q[1]);
  

  //new
  List return_list;

  if(need_mean){
    //posterior mean
    MatrixXd s_1st=Get_s_1st(m_a[0],m_a[1], C_R_K_Q[0],S_KK[1]);
    
    return_list.push_back(s_1st);
  }
  
  MatrixXd S_matrix;
  List S_list=S_KK[0];
  List C=C_R_K_Q[0];
  Eigen::MatrixXd C_here=C[n]; //it looks C has n+1 terms but m only have m terms?
  
  
  if(need_sample){
  
    //posterior sample
    VectorXd sample_theta_1_record=Eigen::VectorXd::Zero(n);
    Eigen::VectorXd theta_sample_here;
    
    //Eigen::VectorXd a;
    
    //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
    Eigen::VectorXd theta_minus_a; 
    
    Eigen::VectorXd random_norm=Eigen::VectorXd::Zero(3);
    for(int t=0;t<3;t++){
      random_norm(t)= R::rnorm(0,1.0);
    }
    
    //the first sample
    
    List m=m_a[0];
    
    Eigen::VectorXd theta_sample;
    
    Eigen::VectorXd m_here=m[n-1];
    
    LLT<MatrixXd> lltOfC_here(C_here);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L_C_here = lltOfC_here.matrixL();   //retrieve factor L  in the decomposition
    
    theta_sample_here=m_here+L_C_here*random_norm;
    sample_theta_1_record[n-1]=theta_sample_here[0];
    
    //done for sample 1
    List KK_list=S_KK[1];
    
    List a=m_a[1];
    VectorXd a_here;
    //VectorXd m_here;
    
    VectorXd h;
    MatrixXd KK_matrix;
    MatrixXd S_matrix_plus_1;
    
    MatrixXd H_matrix;
    /*
     for(int t=n-2;t>=0;t--){
     KK_matrix=KK[t];
     a_vec=a[t+1];
     m_vec=m[t];
     s=m_vec+KK_matrix.transpose()*(s-a_vec);
     s_1st[t]=s[0];
     }
     */
    for(int t=n-2;t>=0;t--){
      
      for(int t_norm=0;t_norm<3;t_norm++){
        random_norm(t_norm)= R::rnorm(0,1.0);
      }
      
      KK_matrix=KK_list[t];
      a_here=a[t+1];
      m_here=m[t];
      h=m_here+KK_matrix.transpose()*(theta_sample_here-a_here);
      
      //not sure why is this
      S_matrix=S_list[t];
      S_matrix_plus_1=S_list[t+1];
      

      H_matrix=S_matrix-KK_matrix.transpose()*S_matrix_plus_1*KK_matrix;
      
      //LLT<MatrixXd> lltOfS(S_matrix);    // compute the cholesky decomposition of R called lltofR
     // MatrixXd L_S = lltOfS.matrixL();   //retrieve factor L  in the decomposition
      
      LLT<MatrixXd> lltOfH(H_matrix);    // compute the cholesky decomposition of R called lltofR
       MatrixXd L_H = lltOfH.matrixL();   //retrieve factor L  in the decomposition
      
      theta_sample_here=h+L_H*random_norm;
      sample_theta_1_record[t]=theta_sample_here[0];
      
      //it seems this
      //sample_theta_1_record[t]=h[0];
      
    }

    
    /*
     List S= S_KK[0];
     MatrixXd S_matrix;
     VectorXd Var=VectorXd::Zero(n);
     for(int t=0;t<n;t++){
     S_matrix=S[t];
     Var[t]=S_matrix(0,0);
     }
     //Var=((Var.array()+VV)*sigma_2_hat).matrix();
     Var=((Var.array()+VV)).matrix();
     */
    
    return_list.push_back(sample_theta_1_record);
  }
  
  //output some some variance, first term of S
  if(need_var){
    VectorXd vec_record=Eigen::VectorXd::Zero(n);
    vec_record[n-1]=C_here(0,0);
    for(int t=n-2;t>=0;t--){
      S_matrix=S_list[t];
      vec_record[t]=S_matrix(0,0);
    }
    return_list.push_back(vec_record);
    
  }
  
  
  //let me also get Q out 
  //return_list.push_back(C_R_K_Q[3]);
  
  
  return return_list;
  
}



//  Nov 25 separate the sample and mean
// [[Rcpp::export]] 
List Kalman_smoother_mean_sample_known_C_R_K_Q_exp(const List& C_R_K_Q, const Eigen::MatrixXd& W0,const List& GG,const List& W,const double VV,
                                               const VectorXi& index_obs,const VectorXd& output,  const bool need_sample, const bool need_mean,const bool need_var){
  
  int n=output.size();
  
  

  
  List m_a=Get_m_a_pred(index_obs, output, GG,C_R_K_Q[2]);
  
  List S_KK=Get_S_KK(index_obs, GG,C_R_K_Q[0],C_R_K_Q[1]);
  
  
  //new
  List return_list;
  
  if(need_mean){
    //posterior mean
    MatrixXd s_1st=Get_s_1st(m_a[0],m_a[1], C_R_K_Q[0],S_KK[1]);
    
    return_list.push_back(s_1st);
  }
  
  MatrixXd S_matrix;
  List S_list=S_KK[0];
  List C=C_R_K_Q[0];
  Eigen::MatrixXd C_here=C[n]; //it looks C has n+1 terms but m only have m terms?
  
  
  if(need_sample){
    
    //posterior sample
    VectorXd sample_theta_1_record=Eigen::VectorXd::Zero(n);
    Eigen::VectorXd theta_sample_here;
    
    //Eigen::VectorXd a;
    
    //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
    Eigen::VectorXd theta_minus_a; 
    
    Eigen::VectorXd random_norm= Eigen::VectorXd::Zero(1);
    random_norm(0)= R::rnorm(0,1.0);
    
    //  Eigen::VectorXd::Zero(3);
    //for(int t=0;t<3;t++){
    //  random_norm(0)= R::rnorm(0,1.0);
    //}
    
    //the first sample
    
    List m=m_a[0];
    
    Eigen::VectorXd theta_sample;
    
    Eigen::VectorXd m_here=m[n-1];
    
    LLT<MatrixXd> lltOfC_here(C_here);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L_C_here = lltOfC_here.matrixL();   //retrieve factor L  in the decomposition
    
    theta_sample_here=m_here+L_C_here*random_norm;
    sample_theta_1_record[n-1]=theta_sample_here[0];
    
    //done for sample 1
    List KK_list=S_KK[1];
    
    List a=m_a[1];
    VectorXd a_here;
    //VectorXd m_here;
    
    VectorXd h;
    MatrixXd KK_matrix;
    MatrixXd S_matrix_plus_1;
    
    MatrixXd H_matrix;

    for(int t=n-2;t>=0;t--){
      
        random_norm(0)= R::rnorm(0,1.0);

      
      KK_matrix=KK_list[t];
      a_here=a[t+1];
      m_here=m[t];
      h=m_here+KK_matrix.transpose()*(theta_sample_here-a_here);
      
      //not sure why is this
      S_matrix=S_list[t];
      S_matrix_plus_1=S_list[t+1];
      
      
      H_matrix=S_matrix-KK_matrix.transpose()*S_matrix_plus_1*KK_matrix;
      
      //LLT<MatrixXd> lltOfS(S_matrix);    // compute the cholesky decomposition of R called lltofR
      // MatrixXd L_S = lltOfS.matrixL();   //retrieve factor L  in the decomposition
      
      LLT<MatrixXd> lltOfH(H_matrix);    // compute the cholesky decomposition of R called lltofR
      MatrixXd L_H = lltOfH.matrixL();   //retrieve factor L  in the decomposition
      
      theta_sample_here=h+L_H*random_norm;
      sample_theta_1_record[t]=theta_sample_here[0];
      
      //it seems this
      //sample_theta_1_record[t]=h[0];
      
    }
    
    
    /*
     List S= S_KK[0];
     MatrixXd S_matrix;
     VectorXd Var=VectorXd::Zero(n);
     for(int t=0;t<n;t++){
     S_matrix=S[t];
     Var[t]=S_matrix(0,0);
     }
     //Var=((Var.array()+VV)*sigma_2_hat).matrix();
     Var=((Var.array()+VV)).matrix();
     */
    
    return_list.push_back(sample_theta_1_record);
  }
  
  //output some some variance, first term of S
  if(need_var){
    VectorXd vec_record=Eigen::VectorXd::Zero(n);
    vec_record[n-1]=C_here(0,0);
    for(int t=n-2;t>=0;t--){
      S_matrix=S_list[t];
      vec_record[t]=S_matrix(0,0);
    }
    return_list.push_back(vec_record);
    
  }
  
  
  //let me also get Q out 
  //return_list.push_back(C_R_K_Q[3]);
  
  
  return return_list;
  
}



// [[Rcpp::export]] 
VectorXd  Sample_matern_5_2(const List& GG,const List&  W,const Eigen::MatrixXd& C0,const double VV,const Eigen::MatrixXd& normal_samples_theta,Eigen::VectorXd& normal_samples_y){ //3* (n+1) times normal damples
  
  int n=GG.size();
  
  // for matern 2.5
  //Eigen::VectorXd random_norm=Eigen::VectorXd::Zero(3);
  
  //for(int t=0;t<3;t++){
  //  random_norm(t)= R::rnorm(0,1.0);
  //}
  
  LLT<MatrixXd> lltOfM(C0);    // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfM.matrixL();   //retrieve factor L  in the decomposition
  
  Eigen::VectorXd sample_cur=L*normal_samples_theta.col(0); // current sample
  
  Eigen::VectorXd sample_record=Eigen::VectorXd::Zero(n);
  
  Eigen::MatrixXd W_mat;
  Eigen::MatrixXd G_mat;
  
  for(int t=0;t<n;t++){
    W_mat=W[t];
    G_mat=GG[t];
    LLT<MatrixXd> lltOfM(W_mat);    
    L = lltOfM.matrixL();
    
    sample_cur=G_mat*sample_cur+L*normal_samples_theta.col(t+1);
    sample_record[t]=sample_cur(0)+sqrt(VV)*normal_samples_y(t);
  }
  
  return sample_record;
  //MatrixXd sample_dlm=Eigen::MatrixXd::Zero(3,3); //
  
  
}

// [[Rcpp::export]] 
VectorXd  Sample_exp(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV,const Eigen::VectorXd normal_samples_theta,Eigen::VectorXd normal_samples_y){ //3* (n+1) times normal damples
  
  int n=GG.size();
  
  // for matern 2.5
  //Eigen::VectorXd random_norm=Eigen::VectorXd::Zero(3);
  
  //for(int t=0;t<3;t++){
  //  random_norm(t)= R::rnorm(0,1.0);
  //}
  
  LLT<MatrixXd> lltOfM(C0);    // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfM.matrixL();   //retrieve factor L  in the decomposition
  
  Eigen::VectorXd sample_cur=L*normal_samples_theta(0); // current sample
  
  Eigen::VectorXd sample_record=Eigen::VectorXd::Zero(n);
  
  Eigen::MatrixXd W_mat;
  Eigen::MatrixXd G_mat;
  
  for(int t=0;t<n;t++){
    W_mat=W[t];
    G_mat=GG[t];
    LLT<MatrixXd> lltOfM(W_mat);    
    L = lltOfM.matrixL();
    
    sample_cur=G_mat*sample_cur+L*normal_samples_theta(t+1);
    sample_record[t]=sample_cur(0)+sqrt(VV)*normal_samples_y(t);
  }
  
  return sample_record;
  //MatrixXd sample_dlm=Eigen::MatrixXd::Zero(3,3); //
  
  
}




////L_times_z
// [[Rcpp::export]] 
VectorXd L_times_z(const List& GG,const List&  W,const Eigen::MatrixXd& W0,const double VV,
                    const Eigen::VectorXd& z){
  
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  
  int n=z.rows();
  

  List    Q_K;
  
  Q_K=Get_Q_K(GG,W,W0,VV);
  
  
  Eigen::VectorXd Q=Q_K[0];
  Eigen::MatrixXd K=Q_K[1];
  
  //return_list.push_back(log_det_R);
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(n);
  
  Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::VectorXd output=Eigen::VectorXd::Zero(n);
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    output[t]=a[0]+z[t]*sqrt_Q[t];
    
    //Y_minus_a_1[t]=(output[t]-a[0]);
    
    // Y_minus_a_1_scaled_vec[t]=(output[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(output[t]-a[0]);
  }
  
  return output;
}

//new code in Feb 2020 to implement MCMC in C++
//Eigs
// [[Rcpp::export]] 
List Eigs(const Eigen::MatrixXd A, const int k){
  
  DenseSymMatProd<double> op(A);
  
  // Construct eigen solver object, requesting the largest three eigenvalues
  SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, k, 2*k);
  // Initialize and compute
  eigs.init();
  //int nconv =
  eigs.compute();
  
  List eigen_list;
  
  eigen_list.push_back(eigs.eigenvalues());
  eigen_list.push_back(eigs.eigenvectors());
  
  return eigen_list;
  
}



//Get_eigs_1
// [[Rcpp::export]] 
List Get_eigs_1(const MatrixXd& R0_1,double log_beta1,String kernel_type,int d){
  double beta1=exp(log_beta1);
  MatrixXd R;
  int n1=R0_1.cols();
  //int ncols=R0_1.cols();
  if(kernel_type=="exp"){
     R=(-(beta1*R0_1).array().pow(1.0)).exp().matrix();
  }else{//else it is Matern 2.5
    const double cnst = sqrt(5.0);
    Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(n1,n1);
    R = cnst*beta1*R0_1;
    R= ((matOnes + R +
            R.array().pow(2.0).matrix()/3.0).cwiseProduct((-R).array().exp().matrix()));
  }
  
  List eigen_list;
  
  if(d<n1){
      //solve eigen problem
    DenseSymMatProd<double> op(R);
    
    // Construct eigen solver object, requesting the largest three eigenvalues
    SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, d, min(2*d,n1));
    // Initialize and compute
    eigs.init();
    eigs.compute();
    //int nconv = eigs.compute();
    
    //first eigen value and second eigen vectors
    eigen_list.push_back(eigs.eigenvalues());
    eigen_list.push_back(eigs.eigenvectors());
  }else{
    EigenSolver<MatrixXd> es(R);
    
    eigen_list.push_back(es.eigenvalues());
    eigen_list.push_back(es.eigenvectors());
    
  }
  return eigen_list;
  
}

//Construct_Cor_KF
// [[Rcpp::export]] 
List Construct_Cor_KF(const VectorXd& log_beta_eta,const VectorXd& delta_x,const VectorXi& index_obs,String kernel_type){
  double lambda;
  double eta;
  
  List Cor_list;
  if(kernel_type=="exp"){
    lambda=1.0*exp(log_beta_eta[0]);
    eta=exp(log_beta_eta[1]);
    
    Cor_list.push_back(Construct_G_exp( delta_x,  lambda));
    Cor_list.push_back(Construct_W0_exp( 1.0,  lambda));
    Cor_list.push_back(Construct_W_exp(1.0, delta_x, lambda, Cor_list[1]) );
  }else{
    lambda=sqrt(5)*exp(log_beta_eta[0]);
    eta=exp(log_beta_eta[1]);
    
    Cor_list.push_back(Construct_G_matern_5_2( delta_x,  lambda));
    Cor_list.push_back(Construct_W0_matern_5_2( 1.0,  lambda));
    Cor_list.push_back(Construct_W_matern_5_2(1.0, delta_x, lambda, Cor_list[1]));
  }
      
    Cor_list.push_back(Get_C_R_K_Q_pred(index_obs, Cor_list[0],Cor_list[2],
                                        Cor_list[1], eta));
      
    return Cor_list;
}


List Construct_Cor_KF_no_C_R_K_Q(const double log_beta,const VectorXd& delta_x,const VectorXi& index_obs,String kernel_type){
  double lambda;
 // double eta;
  
  List Cor_list;
  if(kernel_type=="exp"){
    lambda=1.0*exp(log_beta);
    //eta=exp(log_beta_eta[1]);
    
    Cor_list.push_back(Construct_G_exp( delta_x,  lambda));
    Cor_list.push_back(Construct_W0_exp( 1.0,  lambda));
    Cor_list.push_back(Construct_W_exp(1.0, delta_x, lambda, Cor_list[1]) );
  }else{
    lambda=sqrt(5)*exp(log_beta);
    //eta=exp(log_beta_eta[1]);
    
    Cor_list.push_back(Construct_G_matern_5_2( delta_x,  lambda));
    Cor_list.push_back(Construct_W0_matern_5_2( 1.0,  lambda));
    Cor_list.push_back(Construct_W_matern_5_2(1.0, delta_x, lambda, Cor_list[1]));
  }
  
  //Cor_list.push_back(Get_C_R_K_Q_pred(index_obs, Cor_list[0],Cor_list[2],
  //                                    Cor_list[1], eta));
  
  return Cor_list;
}

  
double log_approx_ref_prior(const VectorXd& param,double nugget, bool nugget_est, const Eigen::VectorXd& CL,const double a,const double b ){
    
    Eigen::VectorXd beta;
    double nu=nugget;
    int param_size=param.size();
    if(!nugget_est){
      beta= param.array().exp().matrix();
    }else{
      beta=param.head(param_size-1).array().exp().matrix(); 
      nu=exp(param[param_size-1]); //nugget
    }
    double t=CL.cwiseProduct(beta).sum()+nu;
    double part_I=-b*t;
    double part_II= a*log(t);
    return part_I+part_II;
}

/*
//Get_beta_lik_S_2_log_det
VectorXd Get_beta_lik_S_2_log_det(const VectorXd& Y_sample_tilde_cur_i_row_t,const VectorXd& log_beta_eta,double tilde_eta,const List& Cov_here,const bool add_log_prior,const VectorXd& prior_par,int n2,double sigma_2_0_cur){
  
  VectorXd return_list=Eigen::VectorXd::Zero(3);
  List C_R_K_Q=Cov_here[3];
  VectorXd L_inv_y_tilde=Get_L_inv_y(Cov_here[0],tilde_eta, C_R_K_Q[3], C_R_K_Q[2],Y_sample_tilde_cur_i_row_t);
 
  double S_2=(L_inv_y_tilde.array()*L_inv_y_tilde.array()).sum()*tilde_eta; //eigen_R1_values[i_row]*eta[i_row]
  VectorXd Q=C_R_K_Q[3];
  double log_det=Q.array().log().sum();
  
  double log_post;
  if(add_log_prior){
    VectorXd CL2=Eigen::VectorXd::Zero(1);
             CL2[0]=prior_par[1];
              
    double log_JR_prior=log_approx_ref_prior(log_beta_eta,0,true,CL2,prior_par[2],prior_par[3]);
    log_post=n2/2.0*(log(tilde_eta))-1/2.0*log_det-S_2/(2.0*sigma_2_0_cur)+log_JR_prior+log_beta_eta.array().sum();
  }else{
    log_post=n2/2.0*(log(tilde_eta))-1/2.0*log_det-S_2/(2.0*sigma_2_0_cur);
  }
  return_list[0]=log_post;
  return_list[1]=S_2;
  return_list[2]=log_det;
  
  return return_list;
    
}
*/

//Get_beta_lik_S_2_log_det
// [[Rcpp::export]] 
VectorXd Get_beta_lik_S_2_log_det(const VectorXd& Y_sample_tilde_cur_i_row_t,const VectorXd& log_beta_eta,double tilde_eta,const List& Cov_here,const List& C_R_K_Q, const bool add_log_prior,const VectorXd& prior_par,int n2,double sigma_2_0_cur){
  
  VectorXd return_list=Eigen::VectorXd::Zero(3);
  //List C_R_K_Q=Cov_here[3];
  VectorXd L_inv_y_tilde=Get_L_inv_y(Cov_here[0],tilde_eta, C_R_K_Q[3], C_R_K_Q[2],Y_sample_tilde_cur_i_row_t);
  
  double S_2=(L_inv_y_tilde.array()*L_inv_y_tilde.array()).sum()*tilde_eta; //eigen_R1_values[i_row]*eta[i_row]
  VectorXd Q=C_R_K_Q[3];
  double log_det=Q.array().log().sum();
  
  double log_post;
  if(add_log_prior){
    VectorXd CL2=Eigen::VectorXd::Zero(1);
    CL2[0]=prior_par[1];
    
    double log_JR_prior=log_approx_ref_prior(log_beta_eta,0,true,CL2,prior_par[2],prior_par[3]);
    log_post=n2/2.0*(log(tilde_eta))-1/2.0*log_det-S_2/(2.0*sigma_2_0_cur)+log_JR_prior+log_beta_eta.array().sum();
  }else{
    log_post=n2/2.0*(log(tilde_eta))-1/2.0*log_det-S_2/(2.0*sigma_2_0_cur);
  }
  return_list[0]=log_post;
  return_list[1]=S_2;
  return_list[2]=log_det;
  
  return return_list;
  
}

//
// [[Rcpp::export]] 
List GOLF_2D_MCMC(const MatrixXd& Y_sample, const MatrixXd & R0_1,const VectorXd&  delta_x_2, String kernel_type,
               const int M, const int M_0, const int d, const VectorXi& missing_index, const MatrixXi& missing_index_all_mat,
               const VectorXd& prior_par,const VectorXd& initial_kernel_par,
               const VectorXd& step_size_kernel_par,const int seed_here,
               const bool need_interval,const double interval_prop, const bool have_mean, const MatrixXd& H2 ){
  
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed_here);  
  
  int n1=Y_sample.rows();
  int n2=Y_sample.cols();
  
  MatrixXd Y_sample_here=Y_sample;
  
  VectorXd log_beta_eta_cur=initial_kernel_par.head(2*d);
  double log_beta1_cur=initial_kernel_par[2*d];
  
  
  List eigen_R1_cur=Get_eigs_1(R0_1,log_beta1_cur,kernel_type,d);
  
  
  
  VectorXd eigen_R1_values_cur=eigen_R1_cur[0];
  MatrixXd eigen_R1_vectors_cur=eigen_R1_cur[1];
  
 // Rcout << "eigen_R1_values_cur" << eigen_R1_values_cur << "\n";
  
 
  VectorXi index_obs= VectorXi::Ones(n2);
  
  VectorXd tilde_eta_vec_cur=(log_beta_eta_cur.segment(d,d).array().exp()/eigen_R1_values_cur.array()).matrix();

  List Cov_2_cur(d);

  VectorXd log_beta_eta_tilde=VectorXd::Zero(2);;
  List Cov_here;
  List C_R_K_Q_cur(d);
  for(int i_d=0;i_d<d;i_d++){
     log_beta_eta_tilde[0]=log_beta_eta_cur[i_d];
     log_beta_eta_tilde[1]=log(tilde_eta_vec_cur[i_d]);
    
     Cov_here=Construct_Cor_KF_no_C_R_K_Q(log_beta_eta_cur[i_d],delta_x_2,index_obs,kernel_type);
    
     Cov_2_cur[i_d]=Cov_here;
     
     C_R_K_Q_cur[i_d]=(Get_C_R_K_Q_pred(index_obs, Cov_here[0],Cov_here[2],
                                        Cov_here[1], tilde_eta_vec_cur[i_d]));
     
  }
    
  //mean 
  //MatrixXd H2=MatrixXd::Ones(n2,1);
  MatrixXd mean_cur=MatrixXd::Zero(n1,n2);
  int q2=H2.cols();
  MatrixXd  B2_cur=MatrixXd::Zero(q2,n1);  
  if(have_mean){
     MatrixXd mean_cur=(H2*B2_cur).transpose();
  }
  double Trace_Y_sample_2=(Y_sample_here-mean_cur).array().pow(2.0).sum();
  MatrixXd Y_sample_tilde_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here-mean_cur);
  double Trace_Y_sample_tilde_2_cur=Y_sample_tilde_cur.array().pow(2.0).sum();
  
  //initialize the post
  double sigma_2_0_cur=0.1;
  VectorXd log_beta_eta=VectorXd::Zero(2);
  MatrixXd beta_post_cur=MatrixXd::Zero(d,3);
  bool add_log_prior=true;
  for(int i_d=0;i_d<d;i_d++){
 //int i_d=0;
    log_beta_eta[0]=log_beta_eta_cur[i_d];
    log_beta_eta[1]=log_beta_eta_cur[i_d+d];
    
    beta_post_cur.row(i_d)=(Get_beta_lik_S_2_log_det(Y_sample_tilde_cur.row(i_d).transpose(),log_beta_eta,tilde_eta_vec_cur[i_d],Cov_2_cur[i_d],C_R_K_Q_cur[i_d],add_log_prior,prior_par,n2,sigma_2_0_cur)).transpose();


  }
  
  
  double tau_sample=R::rgamma((n1)*n2/2.0,2.0/(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+beta_post_cur.col(1).array().sum()) );
    
    sigma_2_0_cur=1.0/tau_sample;
      
      
  //set something to store and uppack some par
  
  //MatrixXd Y_sample_mean_all_record=MatrixXd::Zero(20,n1*n2);
  int index_print=0;
  
  MatrixXd param_record=MatrixXd::Zero(M,2*d+2);

  //variables for the section 1 (sampling beta eta)
  
  double log_beta_propose;
  double log_eta_propose;
  double  tilde_eta_propose;
  //VectorXd log_beta_eta_tilde_propose=VectorXd::Zero(2);
  VectorXd log_beta_eta_propose=VectorXd::Zero(2);
  VectorXd beta_post_S_2_log_det_propose;
  List Cov_2_propose;
  
  double r_ratio;
  bool decision;
  double u;
  
  
  //variables for the section 2 (sampling beta_1)
  double log_JR_prior1_cur;
  double log_post_beta1_cur;
  VectorXd CL1=Eigen::VectorXd::Zero(1);
  CL1[0]=prior_par[0];
  double log_beta1_propose;
  List eigen_R1_propose;
  VectorXd eigen_R1_values_propose;
  MatrixXd eigen_R1_vectors_propose;
  List C_R_K_Q_propose(d);
  List C_R_K_Q_beta1_propose(d);
  VectorXd tilde_eta_vec_propose;
  VectorXd log_beta1_cur_vec=VectorXd::Zero(1);
  MatrixXd Y_sample_tilde_propose; 
  double Trace_Y_sample_tilde_2_propose;
  List C_R_K_Q_here;
  MatrixXd beta_post_propose=MatrixXd::Zero(d,3);
  VectorXd Q_here;
  double log_JR_prior;
  VectorXd L_inv_y_tilde;
  VectorXd CL2=Eigen::VectorXd::Zero(1);
  CL2[0]=prior_par[1];
  VectorXd log_beta1_propose_vec=Eigen::VectorXd::Zero(1);
  double log_JR_prior1_propose;
  double log_post_beta1_propose;
  List C_R_K_Q_beta1_propose_here;
  
  //section 3, sample sigma_2_0
  VectorXd sigma_2_cur;
  double sigma_2_0_cur_sqrt;
  //section 4, sample B2
  //if(have_mean){
    MatrixXd L_inv_H2_var=MatrixXd::Zero(n2,q2);
    MatrixXd A_T_Y_sample_cur;  
    VectorXd L_inv_y_var_2;
    MatrixXd H2_t_R_inv_H2;
    VectorXd H2_t_R_inv_y;
    MatrixXd H2_t_R_inv_H2_inv;
    MatrixXd L_B2_tilde_sample;
    //VectorXd B2_tilde_sample_hat;
    MatrixXd tilde_B2_sample=MatrixXd::Zero(q2,d);
    MatrixXd tilde_B2_hat=MatrixXd::Zero(q2,d);
    //LLT<MatrixXd>  lltOfH2;
    VectorXd normal_sample_q2=VectorXd::Zero(q2);
    MatrixXd normal_sample_q2_n1_mat=MatrixXd::Zero(q2,n1);
    
    MatrixXd H2_t_H2_inv_H2_t_Y;
    MatrixXd B2_hat_c;
    MatrixXd B2_hat;
    
    MatrixXd H2_t_H2_inv=(H2.transpose()*H2).inverse();
      
    MatrixXd H2_t_H2_inv_H2_t=H2_t_H2_inv*H2.transpose();
    LLT<MatrixXd> lltOfH2t_H2_inv(H2_t_H2_inv);
    MatrixXd L_H2_t_H2_inv = lltOfH2t_H2_inv.matrixL();
    MatrixXd part1;
    MatrixXd B2_cur_c;
  //}
  //section 5
  bool need_mean;
  bool need_sample;
  bool need_var;
  MatrixXd z_sample_matrix=MatrixXd::Zero(d,n2);
  MatrixXd Y_sample_mean_all=MatrixXd::Zero(n1,n2);
  MatrixXd Y_sample_all;
  //VectorXd normal_sample_y=VectorXd::Zero(n1,n2);
  List KF_smoother_sample_list;
  VectorXd KF_smoother_sample_vec;
  VectorXd eigen_value_sigma_2_sqrt;
  MatrixXd Y_sample_mean_sum=MatrixXd::Zero(n1,n2);
  
  //std::random_device rd{};
  //std::mt19937 generator{rd()};
  //srand(1);
  
  //std::random_device rand();
  //can change this seed number. Now I set it to be 1
  //std::mt19937 generator{1};
  std::mt19937 generator(seed_here);
  
  //std::default_random_engine generator;
  
  std::normal_distribution<double> distribution_norm{0,1.0};
  std::uniform_real_distribution<double> distribution_unif(0, 1.0);
  
  //define interval
  int n_missing=missing_index_all_mat.rows();
  MatrixXd lower_interval;
  MatrixXd upper_interval;
  MatrixXd initial_interval;
  VectorXd initial_interval_row_vec;
  
  int M_interval;
  VectorXd Y_sample_missing=VectorXd::Zero(n_missing);
  int index_interval;

  if(need_interval){
    M_interval=(int) (M-M_0)*interval_prop;
    lower_interval=MatrixXd::Zero(n_missing,M_interval);
    upper_interval=MatrixXd::Zero(n_missing,M_interval);
    initial_interval=MatrixXd::Zero(n_missing,M_interval);
  }
  
  //let me output all sample mssing, should comment after test
  //MatrixXd Y_sample_missing_record=MatrixXd::Zero(n_missing,M-M_0);
  
  //record eigenvalues
  
  VectorXd eigenvalues_sum=VectorXd::Zero(M);
    
    
  //start of MCMC
  for(int i_M=0; i_M<M;i_M++){
  //  int i_M=0;
  //Rcout << " i_M " << i_M << "log_beta1_cur"<<log_beta1_cur <<"\n";
  
    
    //section 1, sample beta and eta
    
    for(int i_row=0;i_row<d;i_row++){
   //int i_row=0;
      //log_beta_propose=log_beta_eta_cur[i_row]+step_size_kernel_par[0]*R::rnorm(0,1.0);
      //log_eta_propose=log_beta_eta_cur[i_row+d]+step_size_kernel_par[1]*R::rnorm(0,1.0); 
      
      log_beta_propose=log_beta_eta_cur[i_row]+step_size_kernel_par[0]*(distribution_norm(generator));
      log_eta_propose=log_beta_eta_cur[i_row+d]+step_size_kernel_par[1]*(distribution_norm(generator)); 
      
      

      tilde_eta_propose=exp(log_eta_propose)/eigen_R1_values_cur[i_row];
      Cov_2_propose=Construct_Cor_KF_no_C_R_K_Q(log_beta_propose,delta_x_2,index_obs,kernel_type);
      
      C_R_K_Q_propose[i_row]=(Get_C_R_K_Q_pred(index_obs, Cov_2_propose[0],Cov_2_propose[2],
                                        Cov_2_propose[1], tilde_eta_propose));
      
        
      log_beta_eta_propose[0]=log_beta_propose;
      log_beta_eta_propose[1]=log_eta_propose;

      
      beta_post_S_2_log_det_propose=(Get_beta_lik_S_2_log_det(Y_sample_tilde_cur.row(i_row).transpose(),log_beta_eta_propose,tilde_eta_propose,Cov_2_propose,C_R_K_Q_propose[i_row],add_log_prior,prior_par,n2,sigma_2_0_cur)).transpose();
      
      
        
      r_ratio=exp(beta_post_S_2_log_det_propose[0]-beta_post_cur(i_row,0) );
        
      //decision=Accept_proposal(r_ratio);
      if(r_ratio>=1){
        decision=true;
      }else{
         u= distribution_unif(generator); 
        if(u<r_ratio){
          decision=true;
        }else{
          decision=false;
        }
      }
      
      if(decision){
        log_beta_eta_cur[i_row]=log_beta_propose;
        log_beta_eta_cur[i_row+d]=log_eta_propose;
        tilde_eta_vec_cur[i_row]=tilde_eta_propose;
        Cov_2_cur[i_row]=Cov_2_propose;
        C_R_K_Q_cur[i_row]=C_R_K_Q_propose[i_row];
        beta_post_cur.row(i_row)=beta_post_S_2_log_det_propose;
      }
      param_record(i_M,i_row)=exp(  log_beta_eta_cur[i_row]);
      param_record(i_M,i_row+d)=exp( log_beta_eta_cur[i_row+d]);
      
    }
    
    //section 2 sample beta1
    
    log_beta1_cur_vec(0)=log_beta1_cur;
    log_JR_prior1_cur=log_approx_ref_prior(log_beta1_cur_vec,0,false,CL1,prior_par[2],prior_par[3]);
    
    log_post_beta1_cur=n2/2.0*tilde_eta_vec_cur.array().log().sum()-1/2.0*beta_post_cur.col(2).sum()-(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+beta_post_cur.col(1).array().sum())/(2.0*sigma_2_0_cur)+log_JR_prior1_cur+log_beta1_cur;

    //new
    //log_beta1_propose=log_beta1_cur+step_size_kernel_par[2]*R::rnorm(0,1.0); 
    log_beta1_propose=log_beta1_cur+step_size_kernel_par[2]*(distribution_norm(generator)); 
    

    eigen_R1_propose=Get_eigs_1(R0_1,log_beta1_propose,kernel_type,d);
    eigen_R1_values_propose=eigen_R1_propose[0];
    eigen_R1_vectors_propose=eigen_R1_propose[1];
    tilde_eta_vec_propose=(log_beta_eta_cur.segment(d,d).array().exp()/eigen_R1_values_propose.array()).matrix();
    
    Y_sample_tilde_propose=eigen_R1_vectors_propose.transpose()*(Y_sample_here-mean_cur);
    Trace_Y_sample_tilde_2_propose=Y_sample_tilde_propose.array().pow(2.0).sum();
    
    for(int i_row=0;i_row<d;i_row++){
    //int i_row=0;
      Cov_here= Cov_2_cur[i_row];
      C_R_K_Q_beta1_propose[i_row]=Get_C_R_K_Q_pred(index_obs, Cov_here[0],Cov_here[2],
                                              Cov_here[1], tilde_eta_vec_propose[i_row]);
      

      C_R_K_Q_beta1_propose_here= C_R_K_Q_beta1_propose[i_row];
      L_inv_y_tilde=Get_L_inv_y(Cov_here[0],tilde_eta_vec_propose[i_row], C_R_K_Q_beta1_propose_here[3], C_R_K_Q_beta1_propose_here[2],Y_sample_tilde_propose.row(i_row).transpose());
      
      
      beta_post_propose(i_row,1)=(L_inv_y_tilde.array()*L_inv_y_tilde.array()).sum()*tilde_eta_vec_propose[i_row];
      
      Q_here=C_R_K_Q_beta1_propose_here[3];
      beta_post_propose(i_row,2)=Q_here.array().log().sum();
        
      log_beta_eta[0]=log_beta_eta_cur[i_row];
      log_beta_eta[1]=log_beta_eta_cur[i_row+d];
        
      log_JR_prior=log_approx_ref_prior(log_beta_eta,0,true,CL2,prior_par[2],prior_par[3]);
        
      beta_post_propose(i_row,0)=n2/2.0*(log(tilde_eta_vec_propose[i_row]))-1.0/2.0*beta_post_propose(i_row,2)-beta_post_propose(i_row,1)/(2.0*sigma_2_0_cur)+log_JR_prior+log_beta_eta.array().sum();

    }
    
    
    log_beta1_propose_vec(0)=log_beta1_propose;
    
    log_JR_prior1_propose=log_approx_ref_prior(log_beta1_propose_vec,0,false,CL1,prior_par[2],prior_par[3]);
    
    log_post_beta1_propose=n2/2.0*tilde_eta_vec_propose.array().log().sum()-1.0/2.0*beta_post_propose.col(2).sum()-(Trace_Y_sample_2-Trace_Y_sample_tilde_2_propose+beta_post_propose.col(1).array().sum())/(2.0*sigma_2_0_cur)+log_JR_prior1_propose+log_beta1_propose;
    
    
    r_ratio=exp(log_post_beta1_propose-log_post_beta1_cur );
    //decision=Accept_proposal(r_ratio);
    
    if(r_ratio>=1){
      decision=true;
    }else{
      u= distribution_unif(generator); 
      if(u<r_ratio){
        decision=true;
      }else{
        decision=false;
      }
    }
    

    if(decision){
      log_beta1_cur=log_beta1_propose;
      C_R_K_Q_cur=clone(C_R_K_Q_beta1_propose); //this needs to be clone
      eigen_R1_cur=eigen_R1_propose;
      eigen_R1_values_cur=eigen_R1_values_propose;
      eigen_R1_vectors_cur=eigen_R1_vectors_propose;

      tilde_eta_vec_cur=tilde_eta_vec_propose;
        
      beta_post_cur=beta_post_propose;
        
      log_post_beta1_cur=log_post_beta1_propose;
        
      Y_sample_tilde_cur=Y_sample_tilde_propose;
        
      //need to think why these two are  the same
      Trace_Y_sample_tilde_2_cur=Trace_Y_sample_tilde_2_propose;
        
    }
    //record eigen_value
    eigenvalues_sum[i_M]=eigen_R1_values_cur.array().sum();
      
    param_record(i_M,2*d)=exp(log_beta1_cur);
    
    //section 3 sample variance
    
    tau_sample=R::rgamma((n1)*n2/2.0,2.0/(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+beta_post_cur.col(1).array().sum()) );
    
    sigma_2_0_cur=1.0/tau_sample;
    param_record(i_M,2*d+1)=sigma_2_0_cur;
    

    sigma_2_cur=(sigma_2_0_cur/(log_beta_eta_cur.segment(d,d).array().exp())).matrix();
    
    sigma_2_0_cur_sqrt=sqrt(sigma_2_0_cur);
   //section 4, sample B2
   if(have_mean){
      A_T_Y_sample_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here);
       
      for(int i_row=0;i_row<d;i_row++){
       // int i_row=0;
        Cov_here=Cov_2_cur[i_row];
        C_R_K_Q_here=C_R_K_Q_cur[i_row];
        
        for(int j_q=0;j_q<q2;j_q++){
  
          L_inv_H2_var.col(j_q)=Get_L_inv_y(Cov_here[0],tilde_eta_vec_cur[i_row], C_R_K_Q_here[3], C_R_K_Q_here[2],H2.col(j_q));
          //normal_sample_q2[j_q]=R::rnorm(0,1.0);
          normal_sample_q2[j_q]=(distribution_norm(generator));
        }
        
        L_inv_y_var_2=Get_L_inv_y(Cov_here[0],tilde_eta_vec_cur[i_row], C_R_K_Q_here[3], C_R_K_Q_here[2], A_T_Y_sample_cur.row(i_row).transpose());
        
  
        H2_t_R_inv_H2 =L_inv_H2_var.transpose()*L_inv_H2_var;
        H2_t_R_inv_y=L_inv_H2_var.transpose()*L_inv_y_var_2;
        H2_t_R_inv_H2_inv=H2_t_R_inv_H2.inverse();
        
        LLT<MatrixXd> lltOfH2(H2_t_R_inv_H2_inv);
        L_B2_tilde_sample = lltOfH2.matrixL();
        tilde_B2_hat.col(i_row)=H2_t_R_inv_H2_inv*H2_t_R_inv_y;
        tilde_B2_sample.col(i_row)= tilde_B2_hat.col(i_row)+sqrt(sigma_2_cur[i_row]*eigen_R1_values_cur[i_row])*(L_B2_tilde_sample*normal_sample_q2);
          
  
         
      }
    
      for(int i_col=0;i_col<n1;i_col++){
        for(int j_q=0;j_q<q2;j_q++){
          
          //normal_sample_q2_n1_mat(j_q,i_col)=R::rnorm(0,1.0);
          
          normal_sample_q2_n1_mat(j_q,i_col)=(distribution_norm(generator));
          
          
        }
      }
    
      H2_t_H2_inv_H2_t_Y=H2_t_H2_inv_H2_t*Y_sample_here.transpose();
      B2_hat_c=H2_t_H2_inv_H2_t_Y-(H2_t_H2_inv_H2_t_Y*eigen_R1_vectors_cur)*eigen_R1_vectors_cur.transpose();
      B2_hat=tilde_B2_hat.leftCols(d)*eigen_R1_vectors_cur.transpose()+B2_hat_c;
      part1=sigma_2_0_cur_sqrt*L_H2_t_H2_inv*normal_sample_q2_n1_mat;
      B2_cur_c=part1-(part1*eigen_R1_vectors_cur)*eigen_R1_vectors_cur.transpose();
      B2_cur=tilde_B2_sample*eigen_R1_vectors_cur.transpose()+B2_hat_c+B2_cur_c;
      mean_cur=(H2*B2_cur).transpose();
  
      Y_sample_tilde_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here-mean_cur);
   }
   
    //section 5, sample Z and Y_sample
    
    need_mean=false;
    need_sample=true;
    need_var=false;
    

    
    for(int i_row=0;i_row<d;i_row++){     
    //  int i_row=0;
      
      if(kernel_type=="exp"){
        Cov_here=Cov_2_cur[i_row];
         KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_exp(C_R_K_Q_cur[i_row],Cov_here[1],  Cov_here[0],  Cov_here[2],tilde_eta_vec_cur[i_row],
                                                                            index_obs,Y_sample_tilde_cur.row(i_row).transpose()/sqrt(eigen_R1_values_cur[i_row]*sigma_2_cur[i_row]),need_sample,need_mean,need_var);

      }else{ //if not exp then it is matern_5_2
        Cov_here=Cov_2_cur[i_row];
         KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_matern_5_2(C_R_K_Q_cur[i_row],Cov_here[1],  Cov_here[0],  Cov_here[2],tilde_eta_vec_cur[i_row],
                                                                             index_obs,Y_sample_tilde_cur.row(i_row).transpose()/sqrt(eigen_R1_values_cur[i_row]*sigma_2_cur[i_row]),need_sample,need_mean,need_var);
      }
      //z_hat_matrix.row(i_row)=KF_smoother_sample_list[0] 
      
      KF_smoother_sample_vec=KF_smoother_sample_list[0];
      z_sample_matrix.row(i_row)=KF_smoother_sample_vec.transpose(); 
        
    }
    eigen_value_sigma_2_sqrt= (eigen_R1_values_cur.array()*sigma_2_cur.array()).sqrt().matrix();
    for(int i_col=0;i_col<n2;i_col++){
      Y_sample_mean_all.col(i_col)=mean_cur.col(i_col)+eigen_R1_vectors_cur*(eigen_value_sigma_2_sqrt.array()*z_sample_matrix.col(i_col).array()).matrix();
    }
    
    //Y_sample_all=Y_sample_mean_all+sqrt(sigma_2_0_cur )*normal_sample_y;
    
    /*
    if((i_M<M_0)||(!need_interval) ){  
      for(int i_missing=0;i_missing<n_missing; i_missing++){
       // Y_sample_here(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))=Y_sample_all(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1));
       //Y_sample_here(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))=Y_sample_mean_all(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))+sigma_2_0_cur_sqrt*R::rnorm(0,1.0);
        
       Y_sample_here(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))=Y_sample_mean_all(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))+sigma_2_0_cur_sqrt*distribution_norm(generator);
         //normal_sample_y[i_missing];
         
       //Y_sample_mean_all+sqrt(sigma_2_0_cur )*normal_sample_y;
      }
    }else{
      for(int i_missing=0;i_missing<n_missing; i_missing++){
        Y_sample_missing[i_missing]=Y_sample_mean_all(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1));
        Y_sample_here(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))= Y_sample_missing[i_missing]+sigma_2_0_cur_sqrt*distribution_norm(generator);
      }
      
    }
    */
    for(int i_missing=0;i_missing<n_missing; i_missing++){
        Y_sample_missing[i_missing]=Y_sample_mean_all(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))+sigma_2_0_cur_sqrt*distribution_norm(generator);
        Y_sample_here(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))= Y_sample_missing[i_missing];
    }
 

    Y_sample_tilde_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here-mean_cur);

      
      if(i_M>=M_0){
        
        //Y_sample_mean_missing_record.col(i_M-M_0)=Y_sample_mean_missing;  //delete
        
        Y_sample_mean_sum=Y_sample_mean_sum+Y_sample_mean_all;
        
        if((i_M-M_0)<M_interval){
          initial_interval.col(i_M-M_0)=Y_sample_missing;
          //upper_interval.col(i_M-M_0)=Y_sample_missing;
          
          if((i_M-M_0)==(M_interval-1)){
            for(int i_missing=0;i_missing<n_missing; i_missing++){
              //std::sort(lower_interval.row(i_missing).data(), lower_interval.row(i_missing).data()+lower_interval.row(i_missing).size()); // ascending
              initial_interval_row_vec=initial_interval.row(i_missing).transpose();
              std::sort(initial_interval_row_vec.data(), initial_interval_row_vec.data()+initial_interval_row_vec.size()); // ascending
              lower_interval.row(i_missing)=initial_interval_row_vec.transpose();
              upper_interval.row(i_missing)=initial_interval_row_vec.transpose();
              
            }
            // std::sort(lower_interval.data(), lower_interval.data()+lower_interval.size()); // ascending
            // std::sort(upper_interval.data(), upper_interval.data()+upper_interval.size()); // ascending
            
          }
          
        }else{
          
          for(int i_missing=0;i_missing<n_missing; i_missing++){
            //if(i_missing==8 & (i_M-M_0==M_interval) ){
            //Rcout << "The value of the beginning of lower_interval.row(i_missing) : " << lower_interval.row(i_missing) << "\n";
            
            if(Y_sample_missing[i_missing]<lower_interval(i_missing,M_interval-1)){//lower_interval_95[i_missing,M_025] is the largest one
              index_interval=M_interval-2;
              while(Y_sample_missing[i_missing]<lower_interval(i_missing,index_interval) ){
                index_interval=index_interval-1;
                if(index_interval==-1){
                  break;
                }
              }
              
              if(index_interval<(M_interval-2) ){
                VectorXd v_lower=lower_interval.block(i_missing,(index_interval+1),1,(M_interval-index_interval-2)).transpose(); 
                
                lower_interval.block(i_missing,(index_interval+2),1,(M_interval-index_interval-2))= v_lower.transpose();
                
                lower_interval(i_missing,(index_interval+1))=Y_sample_missing[i_missing];
              }else{
                lower_interval(i_missing,M_interval-1)=Y_sample_missing[i_missing];
              }
              
            }
            //}
            
            if(Y_sample_missing[i_missing]>upper_interval(i_missing,0)){//lower_interval_95[i_missing,M_025] is the largest one
              index_interval=1;
              while(Y_sample_missing[i_missing]>upper_interval(i_missing,index_interval) ){
                index_interval=index_interval+1;
                if(index_interval==M_interval){
                  break;
                }
              }
              if(index_interval==1 ){
                upper_interval(i_missing,0)=Y_sample_missing[i_missing];
                
              }else{
                VectorXd v_upper=upper_interval.block(i_missing,1,1,(index_interval-1)).transpose();
                
                upper_interval.block(i_missing,0,1,(index_interval-1))= v_upper.transpose();
                
                //upper_interval(i_missing,0:(index_interval-2) )= upper_interval(i_missing,1:(index_interval-1));
                upper_interval(i_missing,(index_interval-1))=Y_sample_missing[i_missing];
                
              }
              
            }
          }
          
        }
        
      }
      
     Trace_Y_sample_2=(Y_sample_here-mean_cur).array().pow(2.0).sum();
     Trace_Y_sample_tilde_2_cur=Y_sample_tilde_cur.array().pow(2.0).sum();
      
      //section 5 update lik
      
      for(int i_row=0;i_row<d;i_row++){
        log_beta_eta[0]=log_beta_eta_cur[i_row];
        log_beta_eta[1]=log_beta_eta_cur[i_row+d];

        beta_post_cur.row(i_row)=(Get_beta_lik_S_2_log_det(Y_sample_tilde_cur.row(i_row).transpose(),log_beta_eta,tilde_eta_vec_cur[i_row],Cov_2_cur[i_row],C_R_K_Q_cur[i_row],add_log_prior,prior_par,n2,sigma_2_0_cur)).transpose();
        
        
      }
     
  }
  

 
   //return_List

  List return_List(5);
  MatrixXd Y_sample_mean_avg=Y_sample_mean_sum/(M-M_0);
    
  return_List[0]=Y_sample_mean_avg;
  return_List[1]=param_record;
  return_List[2]=eigenvalues_sum;
  
  if(need_interval){
    return_List[3]=lower_interval.col(M_interval-1);
    return_List[4]=upper_interval.col(0);
    //return_List[3]=lower_interval;
    //return_List[4]=upper_interval;
    
    //return_List[4]=Y_sample_missing_record; //should comment out
  }
  

 
  return return_List;
    
    
}


  
//test sort 
/*
VectorXd sort_try(const VectorXd &x){
  VectorXd ans=x;
  std::sort(ans.data(), ans.data()+ans.size());
  return ans;
  
}
*/

// [[Rcpp::export]] 
List GOLF_2D_separable_MCMC(const MatrixXd& Y_sample, const MatrixXd & R0_1,const VectorXd&  delta_x_2, String kernel_type,
                  const int M, const int M_0, const int d, const VectorXi& missing_index, const MatrixXi& missing_index_all_mat,
                  const VectorXd& prior_par,const VectorXd& initial_kernel_par,
                  const VectorXd& step_size_kernel_par,const int seed_here,
                  const bool need_interval,const double interval_prop, const bool have_mean, const MatrixXd& H2 ){
  
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed_here);  
  
  int n1=Y_sample.rows();
  int n2=Y_sample.cols();
  
  MatrixXd Y_sample_here=Y_sample;
  
  VectorXd log_beta_eta_cur=initial_kernel_par.head(2);
  double log_beta1_cur=initial_kernel_par[2];
  
  
  List eigen_R1_cur=Get_eigs_1(R0_1,log_beta1_cur,kernel_type,d);
  
  
  
  VectorXd eigen_R1_values_cur=eigen_R1_cur[0];
  MatrixXd eigen_R1_vectors_cur=eigen_R1_cur[1];
  
  // Rcout << "eigen_R1_values_cur" << eigen_R1_values_cur << "\n";
  
  
  VectorXi index_obs= VectorXi::Ones(n2);
  
  VectorXd tilde_eta_vec_cur=(exp(log_beta_eta_cur[1])/eigen_R1_values_cur.array()).matrix();
  

  List  Cov_2_cur(3);
  Cov_2_cur=Construct_Cor_KF_no_C_R_K_Q(log_beta_eta_cur[0],delta_x_2,index_obs,kernel_type);
  
  List C_R_K_Q_cur(d);
  for(int i_d=0;i_d<d;i_d++){
    C_R_K_Q_cur[i_d]=(Get_C_R_K_Q_pred(index_obs, Cov_2_cur[0],Cov_2_cur[2],
                                       Cov_2_cur[1], tilde_eta_vec_cur[i_d]));
  }
  
  
  //mean 
  //MatrixXd H2=MatrixXd::Ones(n2,1);
  MatrixXd mean_cur=MatrixXd::Zero(n1,n2);
  int q2=H2.cols();
  MatrixXd  B2_cur=MatrixXd::Zero(q2,n1);  
  if(have_mean){
    MatrixXd mean_cur=(H2*B2_cur).transpose();
  }
  double Trace_Y_sample_2=(Y_sample_here-mean_cur).array().pow(2.0).sum();
  MatrixXd Y_sample_tilde_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here-mean_cur);
  double Trace_Y_sample_tilde_2_cur=Y_sample_tilde_cur.array().pow(2.0).sum();
  
  
  //initialize the post
  double sigma_2_0_cur=0.1;
  //VectorXd log_beta_eta=VectorXd::Zero(2);
  MatrixXd beta_post_cur=MatrixXd::Zero(d,3);
  bool add_log_prior=false;
  for(int i_d=0;i_d<d;i_d++){
    //int i_d=0;
    //log_beta_eta[0]=log_beta_eta_cur[i_d];
    //log_beta_eta[1]=log_beta_eta_cur[i_d+d];
    
    //note here it should be Cov_2_cur not Cov_2_cur[i_d]
    beta_post_cur.row(i_d)=(Get_beta_lik_S_2_log_det(Y_sample_tilde_cur.row(i_d).transpose(),log_beta_eta_cur,tilde_eta_vec_cur[i_d],Cov_2_cur,C_R_K_Q_cur[i_d],add_log_prior,prior_par,n2,sigma_2_0_cur)).transpose();
    
  }
  
  
  double tau_sample=R::rgamma((n1)*n2/2.0,2.0/(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+beta_post_cur.col(1).array().sum()) );
  
  sigma_2_0_cur=1.0/tau_sample;
  
  
  //set something to store and uppack some par
  
  //MatrixXd Y_sample_mean_all_record=MatrixXd::Zero(20,n1*n2);
  int index_print=0;
  
  MatrixXd param_record=MatrixXd::Zero(M,4);
  
  //variables for the section 1 (sampling beta eta)
  
  double log_beta_propose;
  double log_eta_propose;
  double  tilde_eta_propose;
  //VectorXd log_beta_eta_tilde_propose=VectorXd::Zero(2);
  VectorXd log_beta_eta_propose=VectorXd::Zero(2);
  VectorXd beta_post_S_2_log_det_propose;
  List Cov_2_propose;
  
  double r_ratio;
  bool decision;
  double u;
  
  //variables for the section 2 (sampling beta_1)
  double log_JR_prior1_cur;
  double log_post_beta1_cur;
  double log_JR_prior_cur;
  double log_JR_prior_propose;
  VectorXd CL1=Eigen::VectorXd::Zero(1);
  CL1[0]=prior_par[0];
  double log_beta1_propose;
  List eigen_R1_propose;
  VectorXd eigen_R1_values_propose;
  MatrixXd eigen_R1_vectors_propose;
  List C_R_K_Q_propose(d);
  //List C_R_K_Q_beta1_propose(d);
  VectorXd tilde_eta_vec_propose;
  VectorXd log_beta1_cur_vec=VectorXd::Zero(1);
  MatrixXd Y_sample_tilde_propose; 
  double Trace_Y_sample_tilde_2_propose;
  List C_R_K_Q_here;
  MatrixXd beta_post_propose=MatrixXd::Zero(d,3);
  VectorXd Q_here;
  double log_JR_prior;
  VectorXd L_inv_y_tilde;
  VectorXd CL2=Eigen::VectorXd::Zero(1);
  CL2[0]=prior_par[1];
  VectorXd log_beta1_propose_vec=Eigen::VectorXd::Zero(1);
  double log_JR_prior1_propose;
  double log_post_beta1_propose;
  //List C_R_K_Q_beta1_propose_here;
  
  //section 3, sample sigma_2_0
  double sigma_2_cur;
  double sigma_2_0_cur_sqrt;
  //section 4, sample B2
  //if(have_mean){
  MatrixXd L_inv_H2_var=MatrixXd::Zero(n2,q2);
  MatrixXd A_T_Y_sample_cur;  
  VectorXd L_inv_y_var_2;
  MatrixXd H2_t_R_inv_H2;
  VectorXd H2_t_R_inv_y;
  MatrixXd H2_t_R_inv_H2_inv;
  MatrixXd L_B2_tilde_sample;
  //VectorXd B2_tilde_sample_hat;
  MatrixXd tilde_B2_sample=MatrixXd::Zero(q2,d);
  MatrixXd tilde_B2_hat=MatrixXd::Zero(q2,d);
  //LLT<MatrixXd>  lltOfH2;
  VectorXd normal_sample_q2=VectorXd::Zero(q2);
  MatrixXd normal_sample_q2_n1_mat=MatrixXd::Zero(q2,n1);
  
  MatrixXd H2_t_H2_inv_H2_t_Y;
  MatrixXd B2_hat_c;
  MatrixXd B2_hat;
  
  MatrixXd H2_t_H2_inv=(H2.transpose()*H2).inverse();
  
  MatrixXd H2_t_H2_inv_H2_t=H2_t_H2_inv*H2.transpose();
  LLT<MatrixXd> lltOfH2t_H2_inv(H2_t_H2_inv);
  MatrixXd L_H2_t_H2_inv = lltOfH2t_H2_inv.matrixL();
  MatrixXd part1;
  MatrixXd B2_cur_c;
  //}
  //section 5
  bool need_mean;
  bool need_sample;
  bool need_var;
  MatrixXd z_sample_matrix=MatrixXd::Zero(d,n2);
  MatrixXd Y_sample_mean_all=MatrixXd::Zero(n1,n2);
  MatrixXd Y_sample_all;
  //VectorXd normal_sample_y=VectorXd::Zero(n1,n2);
  List KF_smoother_sample_list;
  VectorXd KF_smoother_sample_vec;
  VectorXd eigen_value_sigma_2_sqrt;
  MatrixXd Y_sample_mean_sum=MatrixXd::Zero(n1,n2);
  
  //std::random_device rd{};
  //std::mt19937 generator{rd()};
  //srand(1);
  
  //std::random_device rand();
  //can change this seed number. Now I set it to be 1
  //std::mt19937 generator{1};
  std::mt19937 generator(seed_here);
  
  //std::default_random_engine generator;
  
  std::normal_distribution<double> distribution_norm{0,1.0};
  std::uniform_real_distribution<double> distribution_unif(0, 1.0);
  
  //define interval
  int n_missing=missing_index_all_mat.rows();
  MatrixXd lower_interval;
  MatrixXd upper_interval;
  MatrixXd initial_interval;
  VectorXd initial_interval_row_vec;
  
  int M_interval;
  VectorXd Y_sample_missing=VectorXd::Zero(n_missing);
  int index_interval;
  
  if(need_interval){
    M_interval=(int) (M-M_0)*interval_prop;
    lower_interval=MatrixXd::Zero(n_missing,M_interval);
    upper_interval=MatrixXd::Zero(n_missing,M_interval);
    initial_interval=MatrixXd::Zero(n_missing,M_interval);
  }
  
  //let me output all sample mssing, should comment after test
  //MatrixXd Y_sample_missing_record=MatrixXd::Zero(n_missing,M-M_0);
  
  //record eigenvalues
  
  VectorXd eigenvalues_sum=VectorXd::Zero(M);
  

  
  
  //start of MCMC
  for(int i_M=0; i_M<M;i_M++){
    //  int i_M=0;
    //Rcout << " i_M " << i_M << "log_beta1_cur"<<log_beta1_cur <<"\n";
    
    
    //section 1 and 2, sample beta and eta, beta1
    
    log_beta_propose=log_beta_eta_cur[0]+step_size_kernel_par[0]*(distribution_norm(generator));
    log_eta_propose=log_beta_eta_cur[1]+step_size_kernel_par[1]*(distribution_norm(generator)); 
    log_beta1_propose=log_beta1_cur+step_size_kernel_par[2]*(distribution_norm(generator)); 
    
    log_beta_eta_propose[0]=log_beta_propose;
    log_beta_eta_propose[1]=log_eta_propose;
      
      

    Cov_2_propose=Construct_Cor_KF_no_C_R_K_Q(log_beta_propose,delta_x_2,index_obs,kernel_type);
    eigen_R1_propose=Get_eigs_1(R0_1,log_beta1_propose,kernel_type,d);
    

    eigen_R1_values_propose=eigen_R1_propose[0];
    eigen_R1_vectors_propose=eigen_R1_propose[1];
    tilde_eta_vec_propose=(exp(log_beta_eta_propose[1])/eigen_R1_values_propose.array()).matrix();
    
    Y_sample_tilde_propose=eigen_R1_vectors_propose.transpose()*(Y_sample_here-mean_cur);
    Trace_Y_sample_tilde_2_propose=Y_sample_tilde_propose.array().pow(2.0).sum();
    
    for(int i_row=0;i_row<d;i_row++){
      //int i_row=0;
      C_R_K_Q_propose[i_row]=Get_C_R_K_Q_pred(index_obs, Cov_2_propose[0],Cov_2_propose[2],
                                              Cov_2_propose[1], tilde_eta_vec_propose[i_row]);
      
      beta_post_propose.row(i_row)=(Get_beta_lik_S_2_log_det(Y_sample_tilde_propose.row(i_row).transpose(),log_beta_eta_propose,tilde_eta_vec_propose[i_row],Cov_2_propose,C_R_K_Q_propose[i_row],add_log_prior,prior_par,n2,sigma_2_0_cur)).transpose();
      
    }
    
    log_beta1_propose_vec(0)=log_beta1_propose;
    log_beta1_cur_vec(0)=log_beta1_cur;
    
    log_JR_prior1_propose=log_approx_ref_prior(log_beta1_propose_vec,0,false,CL1,prior_par[2],prior_par[3]);
    log_JR_prior1_cur=log_approx_ref_prior(log_beta1_cur_vec,0,false,CL1,prior_par[2],prior_par[3]);
    
    log_JR_prior_propose=log_approx_ref_prior(log_beta_eta_propose,0,true,CL2,prior_par[2],prior_par[3]);
    log_JR_prior_cur=log_approx_ref_prior(log_beta_eta_cur,0,true,CL2,prior_par[2],prior_par[3]);


    log_post_beta1_propose=n2/2.0*tilde_eta_vec_propose.array().log().sum()-1.0/2.0*beta_post_propose.col(2).sum()-(Trace_Y_sample_2-Trace_Y_sample_tilde_2_propose+beta_post_propose.col(1).array().sum())/(2.0*sigma_2_0_cur)
      +log_JR_prior_propose+log_JR_prior1_propose+log_beta_eta_propose.sum()+log_beta1_propose;
    
    log_post_beta1_cur=n2/2.0*tilde_eta_vec_cur.array().log().sum()-1.0/2.0*beta_post_cur.col(2).sum()-(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+beta_post_cur.col(1).array().sum())/(2.0*sigma_2_0_cur)
      +log_JR_prior_cur+log_JR_prior1_cur+log_beta_eta_cur.sum()+log_beta1_cur;
    
    r_ratio=exp(log_post_beta1_propose-log_post_beta1_cur );
    //decision=Accept_proposal(r_ratio);
    
    if(r_ratio>=1){
      decision=true;
    }else{
      u= distribution_unif(generator); 
      if(u<r_ratio){
        decision=true;
      }else{
        decision=false;
      }
    }
    
    
    if(decision){
      log_beta1_cur=log_beta1_propose;
      log_beta_eta_cur=log_beta_eta_propose;
      Cov_2_cur=Cov_2_propose;
      C_R_K_Q_cur=clone(C_R_K_Q_propose); //this needs to be clone
      eigen_R1_cur=eigen_R1_propose;
      eigen_R1_values_cur=eigen_R1_values_propose;
      eigen_R1_vectors_cur=eigen_R1_vectors_propose;
      
      tilde_eta_vec_cur=tilde_eta_vec_propose;
      
      beta_post_cur=beta_post_propose;
      
      log_post_beta1_cur=log_post_beta1_propose;
      
      Y_sample_tilde_cur=Y_sample_tilde_propose;
      
      //need to think why these two are  the same
      Trace_Y_sample_tilde_2_cur=Trace_Y_sample_tilde_2_propose;
      
    }
    //record eigen_value
    eigenvalues_sum[i_M]=eigen_R1_values_cur.array().sum();
    
    param_record(i_M,0)=exp(log_beta_eta_cur[0]);
    param_record(i_M,1)=exp(log_beta_eta_cur[1]);
    
    param_record(i_M,2)=exp(log_beta1_cur);
    
    //section 3 sample variance
    
    tau_sample=R::rgamma((n1)*n2/2.0,2.0/(Trace_Y_sample_2-Trace_Y_sample_tilde_2_cur+beta_post_cur.col(1).array().sum()) );
    
    sigma_2_0_cur=1.0/tau_sample;
    param_record(i_M,3)=sigma_2_0_cur;
    
    
    sigma_2_cur=(sigma_2_0_cur/exp(log_beta_eta_cur[1]));
    
    sigma_2_0_cur_sqrt=sqrt(sigma_2_0_cur);
    //section 4, sample B2
    if(have_mean){
      A_T_Y_sample_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here);
      
      for(int i_row=0;i_row<d;i_row++){
        // int i_row=0;
        //Cov_here=Cov_2_cur[i_row];
        C_R_K_Q_here=C_R_K_Q_cur[i_row];
        
        for(int j_q=0;j_q<q2;j_q++){
          
          L_inv_H2_var.col(j_q)=Get_L_inv_y(Cov_2_cur[0],tilde_eta_vec_cur[i_row], C_R_K_Q_here[3], C_R_K_Q_here[2],H2.col(j_q));
          //normal_sample_q2[j_q]=R::rnorm(0,1.0);
          normal_sample_q2[j_q]=(distribution_norm(generator));
        }
        
        L_inv_y_var_2=Get_L_inv_y(Cov_2_cur[0],tilde_eta_vec_cur[i_row], C_R_K_Q_here[3], C_R_K_Q_here[2], A_T_Y_sample_cur.row(i_row).transpose());
        
        
        H2_t_R_inv_H2 =L_inv_H2_var.transpose()*L_inv_H2_var;
        H2_t_R_inv_y=L_inv_H2_var.transpose()*L_inv_y_var_2;
        H2_t_R_inv_H2_inv=H2_t_R_inv_H2.inverse();
        
        LLT<MatrixXd> lltOfH2(H2_t_R_inv_H2_inv);
        L_B2_tilde_sample = lltOfH2.matrixL();
        tilde_B2_hat.col(i_row)=H2_t_R_inv_H2_inv*H2_t_R_inv_y;
        tilde_B2_sample.col(i_row)= tilde_B2_hat.col(i_row)+sqrt(sigma_2_cur*eigen_R1_values_cur[i_row])*(L_B2_tilde_sample*normal_sample_q2);
        
        
        
      }
      
      for(int i_col=0;i_col<n1;i_col++){
        for(int j_q=0;j_q<q2;j_q++){
          
          //normal_sample_q2_n1_mat(j_q,i_col)=R::rnorm(0,1.0);
          
          normal_sample_q2_n1_mat(j_q,i_col)=(distribution_norm(generator));
          
          
        }
      }
      
      H2_t_H2_inv_H2_t_Y=H2_t_H2_inv_H2_t*Y_sample_here.transpose();
      B2_hat_c=H2_t_H2_inv_H2_t_Y-(H2_t_H2_inv_H2_t_Y*eigen_R1_vectors_cur)*eigen_R1_vectors_cur.transpose();
      B2_hat=tilde_B2_hat.leftCols(d)*eigen_R1_vectors_cur.transpose()+B2_hat_c;
      part1=sigma_2_0_cur_sqrt*L_H2_t_H2_inv*normal_sample_q2_n1_mat;
      B2_cur_c=part1-(part1*eigen_R1_vectors_cur)*eigen_R1_vectors_cur.transpose();
      B2_cur=tilde_B2_sample*eigen_R1_vectors_cur.transpose()+B2_hat_c+B2_cur_c;
      mean_cur=(H2*B2_cur).transpose();
      
      Y_sample_tilde_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here-mean_cur);
    }
    
    //section 5, sample Z and Y_sample

    need_mean=false;
    need_sample=true;
    need_var=false;
    
    
    
    for(int i_row=0;i_row<d;i_row++){     
      //  int i_row=0;
      
      if(kernel_type=="exp"){
        //Cov_here=Cov_2_cur[i_row];
        KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_exp(C_R_K_Q_cur[i_row],Cov_2_cur[1],  Cov_2_cur[0],  Cov_2_cur[2],tilde_eta_vec_cur[i_row],
                                                                              index_obs,Y_sample_tilde_cur.row(i_row).transpose()/sqrt(eigen_R1_values_cur[i_row]*sigma_2_cur),need_sample,need_mean,need_var);
        
      }else{ //if not exp then it is matern_5_2
        //Cov_here=Cov_2_cur[i_row];
        KF_smoother_sample_list=Kalman_smoother_mean_sample_known_C_R_K_Q_matern_5_2(C_R_K_Q_cur[i_row],Cov_2_cur[1],  Cov_2_cur[0],  Cov_2_cur[2],tilde_eta_vec_cur[i_row],
                                                                                     index_obs,Y_sample_tilde_cur.row(i_row).transpose()/sqrt(eigen_R1_values_cur[i_row]*sigma_2_cur),need_sample,need_mean,need_var);
      }
      //z_hat_matrix.row(i_row)=KF_smoother_sample_list[0] 
      
      KF_smoother_sample_vec=KF_smoother_sample_list[0];
      z_sample_matrix.row(i_row)=KF_smoother_sample_vec.transpose(); 
      
    }
    eigen_value_sigma_2_sqrt= (eigen_R1_values_cur.array()*sigma_2_cur).sqrt().matrix();
    for(int i_col=0;i_col<n2;i_col++){
      Y_sample_mean_all.col(i_col)=mean_cur.col(i_col)+eigen_R1_vectors_cur*(eigen_value_sigma_2_sqrt.array()*z_sample_matrix.col(i_col).array()).matrix();
    }
    
    //Y_sample_all=Y_sample_mean_all+sqrt(sigma_2_0_cur )*normal_sample_y;
    
    for(int i_missing=0;i_missing<n_missing; i_missing++){
      Y_sample_missing[i_missing]=Y_sample_mean_all(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))+sigma_2_0_cur_sqrt*distribution_norm(generator);
      Y_sample_here(missing_index_all_mat(i_missing,0),missing_index_all_mat(i_missing,1))= Y_sample_missing[i_missing];
    }
    
    
    Y_sample_tilde_cur=eigen_R1_vectors_cur.transpose()*(Y_sample_here-mean_cur);
    
    
    if(i_M>=M_0){
      
      //Y_sample_mean_missing_record.col(i_M-M_0)=Y_sample_mean_missing;  //delete
      
      Y_sample_mean_sum=Y_sample_mean_sum+Y_sample_mean_all;
      
      if((i_M-M_0)<M_interval){
        initial_interval.col(i_M-M_0)=Y_sample_missing;
        //upper_interval.col(i_M-M_0)=Y_sample_missing;
        
        if((i_M-M_0)==(M_interval-1)){
          for(int i_missing=0;i_missing<n_missing; i_missing++){
            //std::sort(lower_interval.row(i_missing).data(), lower_interval.row(i_missing).data()+lower_interval.row(i_missing).size()); // ascending
            initial_interval_row_vec=initial_interval.row(i_missing).transpose();
            std::sort(initial_interval_row_vec.data(), initial_interval_row_vec.data()+initial_interval_row_vec.size()); // ascending
            lower_interval.row(i_missing)=initial_interval_row_vec.transpose();
            upper_interval.row(i_missing)=initial_interval_row_vec.transpose();
            
          }
          // std::sort(lower_interval.data(), lower_interval.data()+lower_interval.size()); // ascending
          // std::sort(upper_interval.data(), upper_interval.data()+upper_interval.size()); // ascending
          
        }
        
      }else{
        
        for(int i_missing=0;i_missing<n_missing; i_missing++){
          //if(i_missing==8 & (i_M-M_0==M_interval) ){
          //Rcout << "The value of the beginning of lower_interval.row(i_missing) : " << lower_interval.row(i_missing) << "\n";
          
          if(Y_sample_missing[i_missing]<lower_interval(i_missing,M_interval-1)){//lower_interval_95[i_missing,M_025] is the largest one
            index_interval=M_interval-2;
            while(Y_sample_missing[i_missing]<lower_interval(i_missing,index_interval) ){
              index_interval=index_interval-1;
              if(index_interval==-1){
                break;
              }
            }
            
            if(index_interval<(M_interval-2) ){
              VectorXd v_lower=lower_interval.block(i_missing,(index_interval+1),1,(M_interval-index_interval-2)).transpose(); 
              
              lower_interval.block(i_missing,(index_interval+2),1,(M_interval-index_interval-2))= v_lower.transpose();
              
              lower_interval(i_missing,(index_interval+1))=Y_sample_missing[i_missing];
            }else{
              lower_interval(i_missing,M_interval-1)=Y_sample_missing[i_missing];
            }
            
          }
          //}
          
          if(Y_sample_missing[i_missing]>upper_interval(i_missing,0)){//lower_interval_95[i_missing,M_025] is the largest one
            index_interval=1;
            while(Y_sample_missing[i_missing]>upper_interval(i_missing,index_interval) ){
              index_interval=index_interval+1;
              if(index_interval==M_interval){
                break;
              }
            }
            if(index_interval==1 ){
              upper_interval(i_missing,0)=Y_sample_missing[i_missing];
              
            }else{
              VectorXd v_upper=upper_interval.block(i_missing,1,1,(index_interval-1)).transpose();
              
              upper_interval.block(i_missing,0,1,(index_interval-1))= v_upper.transpose();
              
              //upper_interval(i_missing,0:(index_interval-2) )= upper_interval(i_missing,1:(index_interval-1));
              upper_interval(i_missing,(index_interval-1))=Y_sample_missing[i_missing];
              
            }
            
          }
        }
        
      }
      
    }
    
    Trace_Y_sample_2=(Y_sample_here-mean_cur).array().pow(2.0).sum();
    Trace_Y_sample_tilde_2_cur=Y_sample_tilde_cur.array().pow(2.0).sum();
    
    //section 5 update lik
    
    for(int i_row=0;i_row<d;i_row++){
      //log_beta_eta[0]=log_beta_eta_cur[i_row];
      //log_beta_eta[1]=log_beta_eta_cur[i_row+d];
      
      beta_post_cur.row(i_row)=(Get_beta_lik_S_2_log_det(Y_sample_tilde_cur.row(i_row).transpose(),log_beta_eta_cur,tilde_eta_vec_cur[i_row],Cov_2_cur,C_R_K_Q_cur[i_row],add_log_prior,prior_par,n2,sigma_2_0_cur)).transpose();
      
      
    }
    
    
  }
  //return_List
  /*
  List return_List(3);
  return_List[0]=B2_cur;
    return_List[1]= mean_cur;
    return_List[2]=  beta_post_cur;
  */
  List return_List(5);
  MatrixXd Y_sample_mean_avg=Y_sample_mean_sum/(M-M_0);
  
  return_List[0]=Y_sample_mean_avg;
  return_List[1]=param_record;
  return_List[2]=eigenvalues_sum;
  
  if(need_interval){
    return_List[3]=lower_interval.col(M_interval-1);
    return_List[4]=upper_interval.col(0);
    //return_List[3]=lower_interval;
    //return_List[4]=upper_interval;
    
    //return_List[4]=Y_sample_missing_record; //should comment out
  }
  
  return return_List;
  
}
  

