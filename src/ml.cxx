#ifndef GUARD_ML_H
#define GUARD_ML_H

#include "../include/ml.hxx"
#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/vectoroperations.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

#include "TMVA/RModel.hxx"
#include "TMVA/RModelParser_Keras.h"
#include "TMVA/RModelParser_PyTorch.h"
// #include "TMVA/RModelParser_ONNX.hxx"
#include "TMVA/SOFIEHelpers.hxx"
#include "TInterpreter.h"
#include "TSystem.h"


#include <iostream>
#include <filesystem>


std::vector<std::vector<double>> readCSV(const std::string& filename) {

  std::vector<std::vector<double>> data;

  // Open the CSV file
  std::ifstream file(filename);

  // Check if the file is open
    if (!file.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return data; // Return an empty vector
    }

    std::string line;

    // Read the file line by line
    while (std::getline(file, line)) {
      std::vector<double> row;
      std::stringstream lineStream(line);
      std::string cell;

      // skip header row
      if (line.find("mean") != std::string::npos)
	continue;

      // Split the line into cells using a comma as a delimiter
      while (std::getline(lineStream, cell, ',')) {
	row.push_back(std::stod(cell.c_str()));
      }

      // Add the row to the 2D vector
      data.push_back(row);
    }

    // Close the file
    file.close();

    return data;
}



namespace ml {


ROOT::RDF::RNode GaussianTransform(ROOT::RDF::RNode df,
				   const std::string &inputname,
				   const std::string &outputname,
				   const std::string &paramfile,
				   const unsigned position,
				   const std::string &vartype) {
  // read params from file
  Logger::get("GaussianTransform")
    ->debug("reading file {}...", paramfile);
  std::vector<std::vector<double>> params = readCSV(paramfile);

  // shifting to index counting
  unsigned var_idx = position - 1;

  auto transform_int = [params, var_idx](const int input){
    double shifted = -10;
    shifted = (input - params[var_idx][1]) / params[var_idx][2];
    Logger::get("GaussianTransform")
    ->debug("transforming var {} with mean {} and std {} from {} to {}", var_idx + 1, params[var_idx][1], params[var_idx][2], input, shifted);
    return shifted;
  };

  auto transform_float = [params, var_idx](const float input){
    double shifted = -10;
    shifted = (input - params[var_idx][1]) / params[var_idx][2];
    Logger::get("GaussianTransform")
    ->debug("transforming var {} with mean {} and std {} from {} to {}", var_idx + 1, params[var_idx][1], params[var_idx][2], input, shifted);
    return shifted;
  };

  auto transform_double = [params, var_idx](const double input){
    double shifted = -10;
    shifted = (input - params[var_idx][1]) / params[var_idx][2];
    Logger::get("GaussianTransform")
    ->debug("transforming var {} with mean {} and std {} from {} to {}", var_idx + 1, params[var_idx][1], params[var_idx][2], input, shifted);
    return shifted;
  };



  if (vartype.rfind("i", 0) == 0) {
    auto df2 = df.Define(outputname, transform_int, {inputname});
    return df2;
  }
  else if (vartype.rfind("d", 0) == 0) {
    auto df2 = df.Define(outputname, transform_double, {inputname});
    return df2;
  }
  else {
    auto df2 = df.Define(outputname, transform_float, {inputname});
    return df2;
  }


}


namespace sofie {

  void CompileModelForRDF(const std::string & headerModelFile, unsigned int ninputs, unsigned int nslots=0) {

    std::string modelName = headerModelFile.substr(0,headerModelFile.find(".hxx"));
    std::string cmd = std::string("#include \"") + headerModelFile + std::string("\"");
    auto ret = gInterpreter->Declare(cmd.c_str());
    if (!ret)
      throw std::runtime_error("Error compiling : " + cmd);
    std::cout << "compiled : " << cmd << std::endl;

    cmd = "auto sofie_functor_" + modelName + " = TMVA::Experimental::SofieFunctor<" + std::to_string(ninputs) + ",TMVA_SOFIE_" + modelName + "::Session>(" + std::to_string(nslots) + ");";
    ret = gInterpreter->Declare(cmd.c_str());
    if (!ret)
      throw std::runtime_error("Error compiling : " + cmd);
    std::cout << "compiled : " << cmd << std::endl;
    std::cout << "Model is ready to be evaluated" << std::endl;
    return;
  }

  std::string SOFIEGenerator(const std::vector<std::string> &input_vec,
                               TMVA::Experimental::SOFIE::RModel &model,
                               const std::string &modelFilePath) {

    std::string modelFileName = std::filesystem::path(modelFilePath).filename();
    std::string::size_type const p(modelFileName.find_last_of('.'));
    std::string modelName = modelFileName.substr(0, p);

    std::string modelHeaderFile = modelName + std::string(".hxx");

    if (!std::filesystem::exists(modelHeaderFile)) {
  Logger::get("SOFIEGenerator")
    ->debug("generating model code...");
  model.Generate();
  Logger::get("SOFIEGenerator")
    ->debug("dumping model code... {}", modelHeaderFile);
  model.OutputGenerated(modelHeaderFile);
  // Logger::get("SOFIEGenerator")
  //   ->debug("printing model code...");
  // model.PrintGenerated();

  std::string modelWeightFile = modelName + std::string(".dat");

  Logger::get("SOFIEGenerator")
    ->debug("compiling model code...");
  CompileModelForRDF(modelHeaderFile, input_vec.size());
    }
    else {
      Logger::get("SOFIEGenerator")
	->debug("model already compiled, skipping");
    }

  std::string sofie_func_str = "sofie_functor_" + modelName + "(rdfslot_, ";
  sofie_func_str += input_vec[0];
  for (unsigned i = 1; i < input_vec.size(); i++) {
    sofie_func_str += ", " + input_vec[i];
  }
  sofie_func_str += ")";

  return sofie_func_str;

}


ROOT::RDF::RNode KerasEvaluate(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vec,
                               const std::string &outputname,
                               const std::string &modelFilePath) {


  Logger::get("KerasEvaluate")
    ->debug("loading model file {} ...", modelFilePath);
  TMVA::Experimental::SOFIE::RModel model = TMVA::Experimental::SOFIE::PyKeras::Parse(modelFilePath);
  Logger::get("KerasEvaluate")
    ->debug("finished loading model");

  auto eval_func = SOFIEGenerator(input_vec, model, modelFilePath);
  Logger::get("KerasEvaluate")
    ->debug("evaluating model code...");

  auto df2 = df.Define(outputname, eval_func);

  return df2;

}


ROOT::RDF::RNode PyTorchEvaluate(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vec,
                               const std::string &outputname,
                               const std::string &modelFilePath) {


  Logger::get("PyTorchEvaluate")
    ->debug("deriving input shapes...");
  std::vector<size_t> inputTensorShapeSequential{1,input_vec.size()};
  std::vector<std::vector<size_t>> inputShapesSequential{inputTensorShapeSequential};

  Logger::get("PyTorchEvaluate")
    ->debug("loading model file {} ...", modelFilePath);
  TMVA::Experimental::SOFIE::RModel model = TMVA::Experimental::SOFIE::PyTorch::Parse(modelFilePath, inputShapesSequential);
  Logger::get("PyTorchEvaluate")
    ->debug("finished loading model");

  auto eval_func = SOFIEGenerator(input_vec, model, modelFilePath);
  Logger::get("PyTorchEvaluate")
    ->debug("evaluating model code...");

  auto df2 = df.Define(outputname, eval_func);

  return df2;

}


// ROOT::RDF::RNode ONNXEvaluate(ROOT::RDF::RNode df,
//                                const std::vector<std::string> &input_vec,
//                                const std::string &outputname,
//                                const std::string &modelFilePath) {


//   bool verboseParser = false;

//   Logger::get("ONNXEvaluate")
//     ->debug("loading model file {} ...", modelFilePath);
//   TMVA::Experimental::SOFIE::RModelParser_ONNX  parser;
//   TMVA::Experimental::SOFIE::RModel model = parser.Parse(modelFilePath, verboseParser);
//   Logger::get("ONNXEvaluate")
//     ->debug("finished loading model");

//   auto eval_func = SOFIEGenerator(input_vec, model, modelFilePath);
//   Logger::get("ONNXEvaluate")
//     ->debug("evaluating model code...");

//   auto df2 = df.Define(outputname, eval_func);

//   return df2;

// }



} // end namespace sofie
} // end namespace ml
#endif /* GUARD_ML_H */
