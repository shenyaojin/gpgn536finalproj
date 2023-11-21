#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

//
// Created by shenyaojin on 11/21/23. shenyaojin@mines.edu
//

using namespace std;
// Function to convert a text file to a 2D vector of floats
vector<vector<float>> datareader(const string& filePath) {
    vector<vector<float>> result;

    ifstream inputFile(filePath);
    if (!inputFile.is_open()) {
        cerr << "Error: Unable to open file " << filePath << endl;
        return result; // Return an empty vector if the file can't be opened
    }

    string line;
    while (getline(inputFile, line)) {
        stringstream ss(line);
        float value;
        vector<float> row;
        while (ss >> value) {
            row.push_back(value);
        }
        result.push_back(row);
    }

    inputFile.close();
    return result;
}
