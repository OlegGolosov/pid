/** @file   Constants.h
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Some constants and enumerators
*/

#ifndef PidConstants_H
#define PidConstants_H 1

#include <unordered_map>

namespace PidParticles {
enum eParticle {
  kBgPos = 1,
  kBgNeg = -1,
  kEPos = 11,
  kENeg = -11,
  kProton = 2212,
  kAntiProton = -2212,
  kPionPos = 211,
  kPionNeg = -211,
  kKaonPos = 321,
  kKaonNeg = -321,
  kDeutron = 1000010020,
  kHe3 = 1000020030,
  kHe4 = 1000020040
};

const std::unordered_map<int, float> masses =
    {
        {kBgPos, 0.},
        {kBgNeg, 0.},
        {kEPos, 0.000511},
        {kENeg, 0.000511},
        {kProton, 0.938},
        {kAntiProton, 0.938},
        {kPionPos, 0.140},
        {kPionNeg, 0.140},
        {kKaonPos, 0.498},
        {kKaonNeg, 0.498},
        {kDeutron, 1.862},
        {kHe3, 2.793},
        {kHe4, 3.724},
};
}// namespace PidParticles

namespace PidFunction {
enum eNames {
  kA = 0,
  kMean,
  kSigma,
  nParams
};
}

#endif
