/** @file   ParticleFit.h
    @class  Pid::ParticleFit
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to store fit resuls for particle specie
*/

#ifndef PidParticleFit_H
#define PidParticleFit_H 1

#include <array>
#include <vector>

#include "Constants.h"
#include "TF1.h"
#include "TCutG.h"

namespace Pid {

class ParticleFit {
 public:
  /**   Default constructor   **/
  ParticleFit();
  explicit ParticleFit(int type) : particle_type_(type){};
  virtual ~ParticleFit() = default;

  [[nodiscard]] std::vector<double> GetFunctionParams(vector<double> x) const;
  double Eval(vector<double> p, double y);

  void SetFitFunction(const TF1& function) { function_ = function; }
  void SetRangeX(vector<float> min, vector<float> max) { minx_ = min, maxx_ = max; }
  void SetRangeY(float min, float max) { miny_ = min, maxy_ = max; }
  void SetParametrizationFunction(uint ivar, const TF1& func) { parametrization_.at(ivar) = func; }
  void SetParametrization(const std::vector<TF1>& parametrization) { 
    if (parametrization.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter functions is not equal to the number of parameters!\nExiting...\n";
      exit(0);
    }
    parametrization_ = parametrization; 
  }
  
  void SetParVariation(const std::vector<std::vector<double>> parVariation) {
    if (parVariation.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter variations is not equal to the number of parameters!\nExiting...\n";
      exit(0);
    }
    parVariation_ = parVariation; 
  }

  void SetParFitLimits(const std::vector<std::vector<double>> parFitLimits) {
    if (parFitLimits.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter fit limits is not equal to the number of parameters!\nExiting...\n";
      exit(0);
    }
    parFitLimits_ = parFitLimits; 
  }

  [[nodiscard]] const TF1& GetFunction() const { return function_; }
  [[nodiscard]] uint GetNpar() const { return function_.GetNpar(); }
  TF1& GetParametrizationFunction(int ipar) { return parametrization_.at(ipar); }
  [[nodiscard]] std::vector<double> GetParVariation(uint ipar) const {return parVariation_.at(ipar);}
  [[nodiscard]] std::vector<double> GetParFitLimits(uint ipar) const {return parFitLimits_.at(ipar);}

  [[nodiscard]] double GetSigma(float p) const { return parametrization_.at(PidFunction::kSigma).Eval(p); }
  [[nodiscard]] double GetMean(float p) const { return parametrization_.at(PidFunction::kMean).Eval(p); }
  [[nodiscard]] double GetIntegral(float p) const { return parametrization_.at(PidFunction::kA).Eval(p) / sqrt(2 * TMath::Pi() / GetSigma(p)); }

  double GetXmin() const { return minx_; }
  double GetXmax() const { return maxx_; }
  double GetYmin() const { return miny_; }
  double GetYmax() const { return maxy_; }

 private:
  TF1 function_;
  std::vector<TF1> parametrization_{};
  std::vector<std::vector<double>> parVariation_{};
  std::vector<std::vector<double>> parFitLimits_{};

  float minx_{-1.};
  float maxx_{-1.};
  float miny_{-1.};
  float maxy_{-1.};

  int npoints_{-1};
  int particle_type_{-1};

  ClassDef(ParticleFit, 2);
};

}// namespace Pid
#endif// PidParticleFit_H
