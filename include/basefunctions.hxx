#ifndef GUARDBASEFUCTIONS_H
#define GUARDBASEFUCTIONS_H

#include "../include/defaults.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/Logger.hxx"
#include "utility/RooFunctorThreadsafe.hxx"
#include "utility/utility.hxx"
#include <nlohmann/json.hpp>

enum Channel { MT = 0, ET = 1, TT = 2, EM = 3 };

namespace basefunctions {

/**
 * @brief Function to filter events based on their run and luminosity block
 * values
 *
 * @param df the dataframe to filter
 * @param json_path the path to the golden json file containing all valid
 runs
 * ans luminosity blocks
 * @param run th column containing the run value
 * @param luminosity the column containing the luminosity block value
 * @param filtername the name of the filter
 * @return a filtered dataframe
 */
inline ROOT::RDF::RNode JSONFilter(ROOT::RDF::RNode df,
                                   const std::string &json_path,
                                   const std::string &run,
                                   const std::string &luminosity,
                                   const std::string &filtername) {
    std::ifstream i(json_path);
    nlohmann::json golden_json;
    i >> golden_json;
    auto jsonFilterlambda = [golden_json](UInt_t run, UInt_t luminosity) {
        bool matched = false;
        // check if the run exists
        if (golden_json.find(std::to_string(run)) != golden_json.end()) {
            // now loop over all luminosity blocks and check if the event is
            // valid
            for (auto &luminosityrange : golden_json[std::to_string(run)]) {
                if (luminosity >= luminosityrange[0] &&
                    luminosity <= luminosityrange[1]) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                Logger::get("JSONFilter")
                    ->debug("Run {} / luminosity {} not in json file", run,
                            luminosity);
            }
        }
        return matched;
    };
    return df.Filter(jsonFilterlambda, {run, luminosity}, filtername);
}

/// Function to add an input quantity under a different name
///
/// \param df the dataframe to add the quantity to
/// \param inputname name of the existing column
/// \param outputname name of the new column
///
/// \returns a dataframe with the new column

template <typename T>
inline ROOT::RDF::RNode rename(ROOT::RDF::RNode df,
                               const std::string &inputname,
                               const std::string &outputname) {
    return df.Define(outputname, [](const T &q) { return q; }, {inputname});
}

/// Function to writeout a variable from a particle. The particle
/// is identified via the index stored in the input vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the variable value
/// \param position index of the position in the input vector
/// \param vecname name of the column containing the input vector
/// \param column name of the column containing the variable values
///
/// \returns a dataframe with the new column

template <typename T>
ROOT::RDF::RNode getvar(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &position, const std::string &vecname,
                        const std::string &column) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &vec, const ROOT::RVec<T> &col) {
            T out = default_value<T>();

            try {
                const int index = vec.at(position);
                out = col.at(index, default_value<T>());
            } catch (const std::out_of_range &e) {
                Logger::get("getvar")->debug(
                    "Index not found, returning dummy value !");
            }

            return out;
        },
        {vecname, column});
}

/// Function to writeout a variable from a RVec index to a NanoAOD column.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the variable value
/// \param position index of the position in the input vector
/// \param column name of the column containing the variable values
///
/// \returns a dataframe with the new column

template <typename T>
ROOT::RDF::RNode getvar(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &position,
                        const std::string &column) {
    return df.Define(outputname,
                     [](const int &pos, const ROOT::RVec<T> &col) {
                         return col.at(pos, default_value<T>());
                     },
                     {position, column});
}

/// Function to add a new quantity with a defined value
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column
/// \param value the value to be added
///
/// \returns a dataframe with the new column

template <typename T>
inline ROOT::RDF::RNode DefineQuantity(ROOT::RDF::RNode df,
                                       const std::string &outputname,
                                       T const &value) {
    return df.Define(outputname, [value]() { return value; }, {});
}

/// Function to sum all elements of the column with name `quantity` with `ROOT::VecOps::RVec<T>`
/// objects.
///
/// This function is a template implementation, i.e. call SumPerEvent<T> if the column
/// `quantity` contains ROOT::VecOps::RVec<T> objects.
///
/// Elements of the `ROOT::VecOps::RVec<T>`, which should enter the sum, can be selected with
/// index lists from the column `collection_index` as `ROOT::VecOps::RVec<int>` objects
/// per entry.
///
/// Internally, `ROOT::VecOps::Sum` is used to calculate the sum. A custom zero element, which is
/// a second optional argument of `ROOT::VecOps::Sum`, can be passed to this function with setting
/// the parameter `zero`. Its default value is `T(0)`. For instance, when dealing with
/// `ROOT::Math::PtEtaPhiMVector` objects, the `zero` parameter must be set to
/// `ROOT::Math::PtEtaPhiMVector(0., 0., 0., 0)` in order to enable summation with this function.
///
/// \param df Input dataframe
/// \param outputname name of the output column
/// \param quantity column name of the vector variable which is summed per entry
/// \param collection_index column name for index lists of the elements which are going to be summed up
/// \param zero zero element passed as second argument to the `ROOT::VecOps::Sum` function
///
/// \returns a dataframe with the new column
template <typename T>
inline ROOT::RDF::RNode SumPerEvent(ROOT::RDF::RNode df,
                                    const std::string &outputname,
                                    const std::string &quantity,
                                    const std::string &collection_index,
                                    const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity, const ROOT::RVec<int> &collection_index) {
        Logger::get("SumPerEvent")->debug(
                    "sum values {} at indices {}", quantity, collection_index);
        T sum = ROOT::VecOps::Sum(ROOT::VecOps::Take(quantity, collection_index), zero);
        Logger::get("SumPerEvent")->debug(
                    "sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity, collection_index});
}

/// Function to sum all elements of the column with name `quantity` with `ROOT::VecOps::RVec<T>`
/// objects.
///
/// This function is a template implementation, i.e. call SumPerEvent<T> if the column
/// `quantity` contains ROOT::VecOps::RVec<T> objects.
///
/// Internally, `ROOT::VecOps::Sum` is used to calculate the sum. A custom zero element, which is
/// a second optional argument of `ROOT::VecOps::Sum`, can be passed to this function with setting
/// the parameter `zero`. Its default value is `T(0)`.
///
/// \param df Input dataframe
/// \param outputname name of the output column
/// \param quantity column name of the vector variable which is summed per entry
/// \param zero zero element passed as second argument to the `ROOT::VecOps::Sum` function
///
/// \returns a dataframe with the new column
template <typename T>
inline ROOT::RDF::RNode SumPerEvent(ROOT::RDF::RNode df,
                                    const std::string &outputname,
                                    const std::string &quantity,
                                    const T zero = T(0)) {
    auto sum_per_event = [zero](const ROOT::RVec<T> &quantity) {
        Logger::get("SumPerEvent")->debug(
                    "sum values {}", quantity);
        T sum = ROOT::VecOps::Sum(quantity, zero);
        Logger::get("SumPerEvent")->debug(
                    "sum {}", sum);
        return sum;
    };
    return df.Define(outputname, sum_per_event, {quantity});
}

/// This function creates a new column `outputname` with the negatives of type of the values
/// in the column `inputname`.
///
/// Note that this function is implemented as a template, so specify the type `T` of the objects
/// in the input column when calling this function with `Negative<T>(...)`.
///
/// \param df Input dataframe
/// \param outputname name of the output column
/// \param inputname column name of the input column
///
/// \returns a dataframe with the new column
template <typename T>
inline ROOT::RDF::RNode Negative(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &inputname) {
    auto negative = [](const T &input) {
        return -input;
    };
    return df.Define(outputname, negative, {inputname});
}

/// This function filters events, where neither of the input flags is true.
/// This is used to filter events which do not pass an underlying requirement in
/// any systematic variation.
///
/// \param df The input dataframe
/// \param filtername The name of the filter, used in the Dataframe report
/// \param flags Parameter pack of column names that contain the considered
/// flags of type bool
///
/// \returns a filtered dataframe
template <class... Flags>
inline ROOT::RDF::RNode FilterFlagsAny(ROOT::RDF::RNode df,
                                       const std::string &filtername,
                                       const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Filter(
        utility::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return Any(flags); }),
        FlagList, filtername);
}

/// This function defines a flag being true if either of the input flags is
/// true.
///
/// \param df The input dataframe
/// \param outputflag The name of the new column
/// \param flags Parameter pack of column names that contain the considered
/// flags of type bool
///
/// \returns a dataframe containing the new column
template <class... Flags>
inline ROOT::RDF::RNode CombineFlagsAny(ROOT::RDF::RNode df,
                                        const std::string &outputflag,
                                        const Flags &...flags) {
    std::vector<std::string> FlagList;
    utility::appendParameterPackToVector(FlagList, flags...);
    const auto nFlags = sizeof...(Flags);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        utility::PassAsVec<nFlags, bool>(
            [](const ROOT::RVec<bool> &flags) { return Any(flags); }),
        FlagList);
}

/// Function to apply a selection on an integer quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param df The input dataframe
/// \param quantity The quantity the selection refers to
/// \param selection The values of the quantity that are accepted
/// \param filtername The name of the filter, used in the Dataframe report
///
/// \returns a filtered dataframe
template <typename T>
inline ROOT::RDF::RNode FilterIntSelection(ROOT::RDF::RNode df,
                                           const std::string &quantity,
                                           const std::vector<T> &selection,
                                           const std::string &filtername) {
    return df.Filter(
        [selection](const T probe) {
            return std::find(selection.begin(), selection.end(), probe) !=
                   selection.end();
        },
        {quantity}, filtername);
}

/// Function to evaluate a `RooWorkspace` function and put the output into a new
/// dataframe column
///
/// \param[in] df The dataframe, where the new column should be added
/// \param[in] outputname name of the new column
/// \param[in] function A `RooFunctor` pointer, which has to be loaded from a
/// Roo Workspace \param[in] inputs a paramter pack containing all column names
/// needed to be able to evaluate the workspace function
///
/// \returns a dataframe with the newly defined output column
template <class... Inputs>
inline ROOT::RDF::RNode
evaluateWorkspaceFunction(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::shared_ptr<RooFunctorThreadsafe> &function,
                          const Inputs &...inputs) {
    Logger::get("evaluateWorkspaceFunction")
        ->debug("Starting evaluation for {}", outputname);
    auto getValue = [function](const ROOT::RVec<float> &values) {
        Logger::get("evaluateWorkspaceFunction")
            ->debug("Type: {} ", typeid(function).name());
        std::vector<double> argvalues(values.begin(), values.end());
        auto result = function->eval(argvalues.data());
        Logger::get("evaluateWorkspaceFunction")->debug("result {}", result);
        return result;
    };
    std::vector<std::string> InputList;
    utility::appendParameterPackToVector(InputList, inputs...);
    const auto nInputs = sizeof...(Inputs);
    Logger::get("evaluateWorkspaceFunction")->debug("nInputs: {} ", nInputs);
    auto df1 = df.Define(
        outputname, utility::PassAsVec<nInputs, float>(getValue), InputList);
    // change back to ROOT::RDF as soon as fix is available
    return df1;
}
/// Helper function to recursively define columns for each entry of a vector
/// quantity
///
/// \param df Input dataframe
/// \param name name of the vector quantity
/// \param names vector of names for the new columns
/// \param idx index of the current recursion loop, should not be set outside
/// this function
///
/// \returns a lambda function to be used in RDF Define
template <typename T>
inline ROOT::RDF::RNode
UnrollVectorQuantity(ROOT::RDF::RNode df, const std::string &name,
                     const std::vector<std::string> &names,
                     const size_t &idx = 0) {
    if (idx >= names.size()) {
        return df;
    }
    auto df1 = df.Define(
        names.at(idx),
        [idx](const std::vector<T> &quantities) { return quantities.at(idx); },
        {name});
    return UnrollVectorQuantity<T>(df1, name, names, idx + 1);
}
/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMax(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = values < cut;
        return mask;
    };
}

/// Function to apply a maximal filter requirement to an integer quantity.
/// Returns true if the value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMaxInt(const int &cut) {
    return [cut](const ROOT::RVec<int> &values) {
        ROOT::RVec<int> mask = values < cut;
        return mask;
    };
}

/// Function to apply a maximal filter requirement to a quantity.
/// Returns true if the absolute value is smaller than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterAbsMax(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = abs(values) < cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to a quantity.
/// Returns true if the value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMin(const float &cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to an integer quantity.
/// Returns true if the value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterMinInt(const int &cut) {
    // As in ROOT, for min we use >=
    return [cut](const ROOT::RVec<int> &values) {
        ROOT::RVec<int> mask = values >= cut;
        return mask;
    };
}

/// Function to apply a minimal filter requirement to a quantity.
/// Returns true if the absolute value is larger than the given cut value
///
/// \param cut The cut value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterAbsMin(const float &cut) {
    return [cut](const ROOT::RVec<float> &values) {
        ROOT::RVec<int> mask = abs(values) >= cut;
        return mask;
    };
}

/// Function to apply an exact filter requirement to an integer quantity.
/// Returns true if the value is equal to the given value
///
/// \param cut The value of the filter
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterEqualInt(const int &cut) {
    return [cut](const ROOT::RVec<int> &values) {
        ROOT::RVec<int> mask = values == cut;
        return mask;
    };
}

/// Function to combine two RVec Masks by multiplying the two RVec elementwise
///
/// \param mask_1 The first mask
/// \param mask_2 The second mask
///
/// \returns a lambda function which returns the multiplication of the two masks
inline auto MultiplyTwoMasks() {
    return [](const ROOT::RVec<Int_t> &mask_1,
              const ROOT::RVec<Int_t> &mask_2) { return mask_1 * mask_2; };
}

/// Function to filter the Tau ID in NanoAOD. The disciminator output is stored
/// as a bitmask. In order to filter based on a given WP, a bitwise comparison
/// between the given ID and the filtered index is performed. Return true if the
/// working point is passed by a tau.
///
/// \param index The bitmask index to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterID(const int &index) {
    return [index](const ROOT::RVec<UChar_t> &IDs) {
        ROOT::RVec<int> mask;
        for (auto const ID : IDs) {
            if (index > 0)
                mask.push_back(std::min(1, int(ID & 1 << (index - 1))));
            else
                mask.push_back(int(1));
        }
        return mask;
    };
}
/// Function to filter the Jet ID in NanoAOD. The Jet ID has 3 possible values
/// (for UL): 0==fail tight ID and fail tightLepVeto, 2==pass tight ID and fail
/// tightLepVeto, 6==pass tight ID and pass tightLepVeto
///
/// \param index The bitmask index to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterJetID(const int &index) {
    return [index](const ROOT::RVec<Int_t> &IDs) {
        ROOT::RVec<int> mask = IDs >= index;
        Logger::get("FilterJetID")->debug("IDs: {}", IDs);
        Logger::get("FilterJetID")->debug("Filtered mask: {}", mask);
        return mask;
    };
}
/// Function to filter the Jet pileup ID in NanoAOD. This ID is applied on jets
/// below a given pt threshold. The jet pileup ID has 4 possible values:
/// 0==fail, 4==pass loose, 6==pass loose and medium, 7==pass loose, medium and
/// tight
///
/// \param PUindex The index to be used for comparison
/// \param PUptcut The pt threshold to be used for comparison
///
/// \returns a lambda function to be used in RDF Define
inline auto FilterJetPUID(const int &PUindex, const float &PUptcut) {
    return [PUindex, PUptcut](const ROOT::RVec<Int_t> &PUIDs,
                              const ROOT::RVec<float> &jet_pts) {
        ROOT::RVec<int> tmp_mask1 = PUIDs >= PUindex;
        ROOT::RVec<int> tmp_mask2 = jet_pts >= PUptcut;
        ROOT::RVec<int> mask = (tmp_mask1 + tmp_mask2) > 0;
        Logger::get("FilterJetPUID")->debug("PUIDs: {}", PUIDs);
        Logger::get("FilterJetPUID")->debug("PUID mask: {}", tmp_mask1);
        Logger::get("FilterJetPUID")->debug("jpts: {}", jet_pts);
        Logger::get("FilterJetPUID")->debug("jetpt mask: {}", tmp_mask2);
        Logger::get("FilterJetPUID")->debug("PUID_final mask: {}", mask);
        return mask;
    };
}
} // namespace basefunctions

#endif /* GUARDBASEFUNCTIONS_H */
