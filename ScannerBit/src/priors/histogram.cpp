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

#include "gambit/ScannerBit/priors/histogram.hpp"

namespace Gambit
{
  namespace Priors
  {

    /// Constructor
    Histogram::Histogram(const std::vector<std::string>& param, const Options& options)
      : BasePrior(param, 1)
      , myparameter(param_names[0])
      {
      // Only valid for 1D parameter transformation
      if (param.size() != 1)
        {
          scan_err << "Invalid input to Histogram prior (in constructor): " << endl
                   << "This prior only works with one-dimensional paramters. "
                   << "Input parameters must be a vector of size 1! (has size " << param.size()
                   << ")" << scan_end;
        }

        // Read the entries we need from the options
        if ( options.hasKey("hist_data") )
        {
          /// Get the histogram from the YAML file
          histogram = options.getValue<std::vector<double>>("hist_data");
          if ( options.hasKey("hist_file") ) {
            scan_err << "You specified both 'hist_data' and 'hist_file' as inputs for the Histogram prior. "
                     << "Please only specify one tp allow for a unique choice." << endl;
          }
        }
        else if ( options.hasKey("hist_file") )
        {
          /// Read column of a data file to get the histogram
          std::string file = options.getValue<std::string>("hist_file");
          ASCIItableReader data (file);
          if (data.getnrow() < 2) {
            scan_err << "The histogram file for the Histogram prior contains only one row. Recall that the values for "
                     << "the histogram need to be provided in a single columns. If you did this intentionally, "
                     << "use the 'fixed_prior' for trivial parameters instead." << endl;
          }
          histogram = data[0];
        }
        else
        {
          scan_err << "You need to specify either the 'hist_data' option or the 'hist_file' option "
                   << "in order for the Histogram prior to work." << endl;
        }

        /// Sort the histogram, get unique entries, and calculate their multiplicies
        std::sort(histogram.begin(), histogram.end());
        unique_entries = histogram;
        n_entries = histogram.size();
        unique_entries.erase(std::unique(unique_entries.begin(), unique_entries.end()), unique_entries.end());
        pmf_values.clear();
        for (auto it = unique_entries.begin(); it != unique_entries.end(); it++)
        {
          auto range = std::equal_range(histogram.begin(), histogram.end(), *it);
          int dist = std::distance(range.first, range.second);
          double pmf = dist / double(n_entries);
          if (pmf <= 0)
          {
            scan_err << "The constructor of the Histogram prior failed to find one of the known histogram "
                     << "entries. This should never happen and is a critical bug; please report it."
                     << endl;
          }
          pmf_values.push_back(pmf);
        }
      }


      /// Transformation from unit interval to the histogram values (inverse CDF)
      void Histogram::transform(const std::vector <double> &unitpars, std::unordered_map <std::string, double> &output) const
      {
        // Only valid for 1D parameter transformation
        if (unitpars.size() != 1)
        {
          scan_err << "Invalid input to Histogram prior (in 'transform'): Input parameters must be a vector of size 1! (has size " << unitpars.size() << ")." << scan_end;
        }

        int index = int(n_entries*unitpars[0]);
        output[myparameter] = histogram[index];
      }


      /// Probability mass function
      double Histogram::operator()(const std::vector<double> &vec) const
      {
        auto entry = std::lower_bound(unique_entries.begin(), unique_entries.end(), vec[0]);
        int index = std::distance(unique_entries.begin(), entry);
        return pmf_values[index];
      }

   }
}
