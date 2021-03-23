//  GAMBIT: Global and Modular BSM Inference Tool
//  *********************************************
///  \file
///
///  Prior function that implements a discrete
///  distribution from a user-provided histogram.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///  \date 2021 March
///
///  *********************************************

#ifndef PRIOR_HISTOGRAM_HPP
#define PRIOR_HISTOGRAM_HPP

#include "gambit/ScannerBit/priors.hpp"
#include "gambit/Utils/yaml_options.hpp"
#include "gambit/Utils/ascii_table_reader.hpp"

#include <vector>
#include <algorithm>

namespace Gambit
{
  namespace Priors
  {
    /// histogram prior for a discrete set of values.
    /// Takes the arguments: [hist_data : hist_file]
    class Histogram : public BasePrior
    {
    private:
      /// Name of the parameter that this prior is supposed to transform
      const std::string &myparameter;

      /// Array that contains the histogram
      std::vector<double> histogram;
      std::vector<double> unique_entries;
      std::vector<double> pmf_values;
      int n_entries;

    public:
      Histogram(const std::vector<std::string>& param, const Options&);

      /// Transformation from unit interval to the histogram values (inverse CDF)
      void transform(const std::vector <double> &unitpars, std::unordered_map <std::string, double> &output) const;
      /// Probability mass function
      double operator()(const std::vector<double> &vec) const;
    };

    LOAD_PRIOR(histogram, Histogram)
  }
}

#endif
