#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <utility>
#include "CN.h"
void writeToCSV(const std::string& filename, const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFile> <outputFile>" << std::endl;
        return 1;
    }

    std::cout << "Started" << std::endl;
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    std::cout << "..." << std::endl;

    std::ifstream input(inputFile);
    if (!input.is_open()) {
        std::cerr << "Error: Unable to open input file " << inputFile << std::endl;
        return 1;
    }

    std::string line;
    std::getline(input, line);  // Read the first line (header)

    std::replace(line.begin(), line.end(), ',', ' ');  // Replace commas with spaces for parsing
    std::istringstream headerStream(line);
    std::string type, option;
    double S_0, K, sigma, Tf, T0;
    int Nx, Nt;

    headerStream >> Tf >> T0 >> S_0 >> sigma >> Nt >> Nx >> K >> option >> type;
    
    // Read the remaining lines for risk-free rates
    std::vector<std::pair<double, double>> r;
    while (std::getline(input, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');  // Replace commas with spaces
        std::istringstream iss(line);
        double first, second;
        if (iss >> first >> second) {
            r.emplace_back(first, second);
        }
        else {
            std::cerr << "Warning: Malformed or empty line in input: " << line << std::endl;
        }
    }
    input.close();
    
    std::cout << "..." << std::endl;
    CN optionPricing(type, option, S_0, K, r, sigma, Nx, Nt, Tf, T0);
    CN americanOption(std::string("American"), option, S_0, K, r, sigma, Nx, Nt, Tf, T0);
    std::cout << "..." << std::endl;
    std::vector<std::vector<double>> optionPrices = optionPricing.solve();
    std::vector<std::vector<double>> delta = optionPricing.delta(0.5);
    std::vector<std::vector<double>> gamma = optionPricing.gamma(0.5);
    std::cout << "..." << std::endl;
    std::vector<std::vector<double>> theta = optionPricing.theta();
    std::vector<std::vector<double>> rho = optionPricing.rho(0.01);
    std::vector<std::vector<double>> vega = optionPricing.vega(0.1);
    std::vector<std::vector<double>> muricanPrices = americanOption.solve();
    std::cout << "..." << std::endl;
    writeToCSV("option_prices.csv", optionPrices);
    writeToCSV("delta.csv", delta);
    writeToCSV("gamma.csv", gamma);
    writeToCSV("theta.csv", theta);
    writeToCSV("rho.csv", rho);
    writeToCSV("vega.csv", vega);
    writeToCSV("USop.csv", muricanPrices);
    std::cout << "..." << std::endl;
    std::string Valid=optionPricing.validate();
    std::string Compare=optionPricing.compare_eu_us_nc();
    std::cout << "..." << std::endl;
    std::ofstream output(outputFile);
    output << Compare;
    output << " , ";
    output << Valid;
    output.close();

    std::cout << "Option pricing and Greeks calculation completed. Results saved to CSV files." << std::endl;
    return 0;
}