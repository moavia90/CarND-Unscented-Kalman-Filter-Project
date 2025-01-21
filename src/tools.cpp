#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
                                 
   VectorXd rmse(4);
	rmse << 0,0,0,0;
   VectorXd res(4);

   for(int i=0; i < estimations.size(); ++i){
      res = estimations[i] - ground_truth[i];
      res = res.array() * res.array();
      rmse += res;
   }

   rmse /=  estimations.size();
   rmse = rmse.array().sqrt();
   return rmse;
}