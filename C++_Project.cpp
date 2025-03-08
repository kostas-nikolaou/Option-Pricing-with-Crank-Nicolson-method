#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <vector>
#include <utility>
#include "CN.h"

void writeToCSV(const std::string& baseFilename, const std::string& suffix, const std::vector<std::vector<double>>& data) {
    std::string filename = baseFilename +"\\" +suffix + ".csv";
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return;
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) file << ",";
        }
        file << "\n";
    }
    file.close();
    std::cout << "Saved: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFile> <outputBasePath>" << std::endl;
        return 1;
    }

    std::cout << "Started" << std::endl;
    std::string inputFile = argv[1];
    std::string outputBasePath = argv[2]; // This will be used as a base for output files

    std::ifstream input(inputFile);
    if (!input.is_open()) {
        std::cerr << "Error: Unable to open input file " << inputFile << std::endl;
        return 1;
    }

    std::string line;
    std::getline(input, line);  // Read the first line (header)

    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream headerStream(line);
    std::string type, option;
    double S_0, K, sigma, Tf, T0, dP, dS, dr;
    int Nx, Nt;

    headerStream >> dr >> dS >> dP >> Tf >> T0 >> S_0 >> sigma >> Nx >> Nt >> K >> type >> option;

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

    CN optionPricing(option, type, S_0, K, r, sigma, Nx, Nt, Tf, T0);
    CN americanOption("American", type, S_0, K, r, sigma, Nx, Nt, Tf, T0);

    std::cout << "Computing values..." << std::endl;

    std::vector<std::vector<double>> optionPrices = optionPricing.solve();
    std::vector<std::vector<double>> delta = optionPricing.delta(dP);
    std::vector<std::vector<double>> gamma = optionPricing.gamma(dP);
    std::vector<std::vector<double>> theta = optionPricing.theta();
    std::vector<std::vector<double>> rho = optionPricing.rho(dr);
    std::vector<std::vector<double>> vega = optionPricing.vega(dS);
    std::vector<std::vector<double>> muricanPrices = americanOption.solve();

    // Writing outputs
    writeToCSV(outputBasePath, "option_prices", optionPrices);
    writeToCSV(outputBasePath, "delta", delta);
    writeToCSV(outputBasePath, "gamma", gamma);
    writeToCSV(outputBasePath, "theta", theta);
    writeToCSV(outputBasePath, "rho", rho);
    writeToCSV(outputBasePath, "vega", vega);
    writeToCSV(outputBasePath, "USop", muricanPrices);

    std::cout << "Final computations..." << std::endl;
    //std::pair<double, int> Valid = optionPricing.validate();
    //double Compare = optionPricing.compare_eu_us_nc();
    std::pair<double, int> Valid = optionPricing.validate();
    double Compare = optionPricing.compare_eu_us_nc();
    std::cout << "..." << std::endl;
    std::ofstream output(outputBasePath + "\\" + "output.csv");
    output << Compare;
    output << ",";
    output << Valid.first;
    output << Valid.second;
    output.close();

    // Writing summary results to a single CSV file
    /*std::ofstream summaryFile(outputBasePath + "\\" + "output.csv");
    if (summaryFile.is_open()) {
        summaryFile << "Compare,Valid.first,Valid.second\n";
        summaryFile << Compare << "," << Valid.first << "," << Valid.second << "\n";
        summaryFile.close();
        std::cout << "Saved: " << outputBasePath << "\\"<<"output.csv" << std::endl;
    }
    else {
        std::cerr << "Error: Unable to open summary file!" << std::endl;
    }
    */
    std::cout << "Option pricing and Greeks calculation completed. Results saved to CSV files." << std::endl;
    return 0;
}
