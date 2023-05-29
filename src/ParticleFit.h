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

using std::string, std::vector, std::cout, std::endl;

class ParticleFit {
 public:
  /**   Default constructor   **/
  ParticleFit();
  explicit ParticleFit(int type) : particle_type_(type){};
  virtual ~ParticleFit() = default;

  [[nodiscard]] std::vector<double> GetInputFunctionParams(float x) const;
  [[nodiscard]] std::vector<double> GetOutputFunctionParams(float x) const;
  double Eval(double x, double y);
  double Integral(double x);

  void SetFitFunction(const TF1& function) { function_ = function; }
  void SetRangeX(float min, float max) { minx_ = min, maxx_ = max; }
  void SetRangeY(float min, float max) { miny_ = min, maxy_ = max; }
  void SetInputParametrizationFunction(uint ivar, const TF1& func) { inputParametrization_.at(ivar) = func; }
  void SetOutputParametrizationFunction(uint ivar, const TF1& func) { outputParametrization_.at(ivar) = func; }
  void SetInputParametrization(const std::vector<TF1>& parametrization) { 
    if (parametrization.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter functions is not equal to the number of input parameters!\nExiting...\n";
      exit(0);
    }
    inputParametrization_ = parametrization; 
  }
  
  void SetOutputParametrization(const std::vector<TF1>& parametrization) { 
    if (parametrization.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter functions is not equal to the number of output parameters!\nExiting...\n";
      exit(0);
    }
    outputParametrization_ = parametrization; 
  }

  void SetParVariation(const std::vector<std::vector<double>> parVariation) {
    if (parVariation.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter variations is not equal to the number of output parameters!\nExiting...\n";
      exit(0);
    }
    parVariation_ = parVariation; 
  }

  void SetParFitOptions(const std::vector<string> parFitOptions)
  {
    if (parFitOptions.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter fit options is not equal to the number of output parameters!\nExiting...\n";
      exit(0);
    }
    parFitOptions_=parFitOptions;
  }

  void SetParFitLimits(const std::vector<std::vector<double>> parFitLimits) 
  {
    if (parFitLimits.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter fit limits is not equal to the number of output parameters!\nExiting...\n";
      exit(0);
    }
    parFitLimits_ = parFitLimits; 
  }

  void SetParLimits(const std::vector<std::vector<double>> parLimits) 
  {
    if (parLimits.size() != GetNpar()) {
      std::cout << "\n\nNumber of parameter limits is not equal to the number of output parameters!\nExiting...\n";
      exit(0);
    }
    parLimits_ = parLimits; 
  }

  void SetType(int type){ particle_type_=type; }

  [[nodiscard]] const TF1& GetFunction() const { return function_; }
  [[nodiscard]] uint GetNpar() const { return function_.GetNpar(); }
  TF1& GetInputParametrizationFunction(int ipar) { return inputParametrization_.at(ipar); }
  TF1& GetOutputParametrizationFunction(int ipar) { return outputParametrization_.at(ipar); }
  vector<TF1> GetOutputParametrization() { return outputParametrization_; }
  [[nodiscard]] std::vector<double> GetParVariation(uint ipar) const {return parVariation_.at(ipar);}
  [[nodiscard]] std::vector<double> GetParLimits(uint ipar) const {return parLimits_.at(ipar);}
  [[nodiscard]] std::vector<double> GetParFitLimits(uint ipar) const {return parFitLimits_.at(ipar);}
  [[nodiscard]] string GetParFitOption(uint ipar) const {return parFitOptions_.at(ipar);}

  [[nodiscard]] double GetSigma(float p) const { return inputParametrization_.at(PidFunction::kSigma).Eval(p); }
  [[nodiscard]] double GetMean(float p) const { return inputParametrization_.at(PidFunction::kMean).Eval(p); }
  [[nodiscard]] double GetIntegral(float p) const { return inputParametrization_.at(PidFunction::kA).Eval(p) / sqrt(2 * TMath::Pi() / GetSigma(p)); }

  [[nodiscard]] double GetXmin() const { return minx_; }
  [[nodiscard]] double GetXmax() const { return maxx_; }
  [[nodiscard]] double GetYmin() const { return miny_; }
  [[nodiscard]] double GetYmax() const { return maxy_; }
  [[nodiscard]] int GetType() const { return particle_type_; }

 private:
  TF1 function_;
  std::vector<TF1> inputParametrization_{};
  std::vector<TF1> outputParametrization_{};
  std::vector<std::vector<double>> parVariation_{};
  std::vector<std::vector<double>> parFitLimits_{};
  std::vector<std::vector<double>> parLimits_{};

  float minx_{-1.};
  float maxx_{-1.};
  float miny_{-1.};
  float maxy_{-1.};

  int npoints_{-1};
  int particle_type_{-1};
  vector<string> parFitOptions_;

  ClassDef(ParticleFit, 2);
};

}// namespace Pid
#endif// PidParticleFit_H
