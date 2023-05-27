#include "ParticleFit.h"

ClassImp(Pid::ParticleFit)

    namespace Pid {

  /// Default constructor   -------------------------------------------
  ParticleFit::ParticleFit() = default;

  /**
* Getter for vector of parameters for a ParticleFit
* @param p track momentum
* @return vector of parameters
*/
  std::vector<double> ParticleFit::GetInputFunctionParams(float x) const {
    std::vector<double> params;

    for (uint i = 0; i < GetNpar(); ++i) {
      params.push_back(inputParametrization_.at(i).Eval(x));
    }

    return params;
  }

  std::vector<double> ParticleFit::GetOutputFunctionParams(float x) const {
    std::vector<double> params;

    for (uint i = 0; i < GetNpar(); ++i) {
      params.push_back(outputParametrization_.at(i).Eval(x));
    }

    return params;
  }

  /**
* Similar to TF1::Eval
* @param x track momentum
* @param y track mass square
* @return vector of parameters
*/
  double ParticleFit::Eval(double x, double y) {
    if (x > maxx_ || x < minx_) return 0.;

    const uint npar = function_.GetNpar();
    if (outputParametrization_.size() != npar)
      exit(1);

    function_.SetParameters(&(GetOutputFunctionParams(x)[0]));
    const double ret = function_.Eval(y);

    return ret;
  }
}
