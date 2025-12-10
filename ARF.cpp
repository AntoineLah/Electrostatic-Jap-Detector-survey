#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include <cmath>

namespace fs = std::filesystem;

struct ROIResult {
    double sum;
    double rate;
    double sumError;
    double rateError;
};

// --- File validation functions ---
bool isValidBackupFile(const std::string& filename) {
    size_t bkupPos = filename.find("bkup");
    if (bkupPos == std::string::npos) return false;
    size_t dotPos = filename.find(".iri", bkupPos);
    if (dotPos == std::string::npos) return false;
    std::string numberStr = filename.substr(bkupPos + 4, dotPos - (bkupPos + 4));
    if (numberStr.empty()) return false;
    for (char c : numberStr) if (!isdigit(c)) return false;
    return true;
}

int extractBackupNumber(const std::string& filename) {
    size_t bkupPos = filename.find("bkup");
    if (bkupPos == std::string::npos) return -1;
    size_t dotPos = filename.find(".iri", bkupPos);
    if (dotPos == std::string::npos) return -1;
    try {
        return std::stoi(filename.substr(bkupPos + 4, dotPos - (bkupPos + 4)));
    } catch (...) {
        return -1;
    }
}

std::string extractParentFolder(const std::string& filepath) {
    fs::path path(filepath);
    return path.parent_path().filename().string();
}

time_t convertToTimestamp(const std::string& foldername) {
    if (foldername.length() != 13 || foldername[6] != '_') return -1;
    int year = 2000 + std::stoi(foldername.substr(0,2));
    int month = std::stoi(foldername.substr(2,2));
    int day = std::stoi(foldername.substr(4,2));
    int hour = std::stoi(foldername.substr(7,2));
    int minute = std::stoi(foldername.substr(9,2));
    int second = std::stoi(foldername.substr(11,2));
    struct tm tm = {0};
    tm.tm_year = year - 1900;
    tm.tm_mon = month - 1;
    tm.tm_mday = day;
    tm.tm_hour = hour;
    tm.tm_min = minute;
    tm.tm_sec = second;
    tm.tm_isdst = -1;
    return mktime(&tm);
}

double readAcquisitionTime(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir " << filepath << std::endl;
        return -1.0;
    }
    double time;
    file >> time;
    file.close();
    return time;
}

std::vector<double> readSpectrum(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir " << filepath << std::endl;
        return {};
    }
    double time;
    file >> time; // ignore acquisition time
    std::vector<double> spectrum(1024,0);
    int bin;
    double counts;
    while(file >> bin >> counts) if(bin < 1024) spectrum[bin] = counts;
    file.close();
    return spectrum;
}

// --- Integrate ROI function ---
std::map<std::string, ROIResult> integrateROIs(const std::vector<double>& spectrum,
                                              const std::vector<std::pair<int,int>>& rois,
                                              double acquisitionTime) {
    std::map<std::string, ROIResult> results;
    for (size_t i=0;i<rois.size();++i) {
        int start = rois[i].first;
        int end = rois[i].second;
        double sum = 0;
        for(int j=start;j<=end;++j) if(j<(int)spectrum.size()) sum += spectrum[j];
        double sumErr = (sum>0) ? std::sqrt(sum) : 0;
        double rate = 0, rateErr=0;
        if(acquisitionTime>0) {
            double factor = 3600.0/acquisitionTime;
            rate = sum*factor;
            rateErr = sumErr*factor;
        }
        results["ROI_"+std::to_string(i+1)] = {sum,rate,sumErr,rateErr};
    }
    return results;
}

// --- Binned TGraphErrors including statistical ---
TGraphErrors* makeBinnedGraph(const std::vector<double>& times,
                              const std::vector<double>& values,
                              int hoursBin,
                              long long tMin, long long tMax) {
    long long binSize = hoursBin*3600;
    int nBins = (tMax-tMin)/binSize + 1;
    std::vector<double> binSums(nBins,0);
    std::vector<int> binCounts(nBins,0);
    std::vector<std::vector<double>> binValues(nBins);

    for(size_t i=0;i<times.size();++i) {
        if(times[i]<tMin || times[i]>tMax) continue;
        int idx = (times[i]-tMin)/binSize;
        binSums[idx] += values[i];
        binCounts[idx] += 1;
        binValues[idx].push_back(values[i]);
    }

    std::vector<double> binTimes, binMeans, binErrors, xErr;
    for(int i=0;i<nBins;++i) {
        if(binCounts[i]==0) continue;
        double mean = binSums[i]/binCounts[i];

        // Statistical error (std dev of bin)
        double statErr = 0;
        if(binCounts[i]>1) {
            double sumSq=0;
            for(double v : binValues[i]) sumSq += (v-mean)*(v-mean);
            statErr = std::sqrt(sumSq)/binCounts[i];
        }

        
        double totalErr = std::sqrt(statErr*statErr);

        binTimes.push_back(tMin + i*binSize + binSize/2.0);
        binMeans.push_back(mean);
        binErrors.push_back(totalErr);
        xErr.push_back((hoursBin*3600/2));
    }

    TGraphErrors* g = new TGraphErrors(binTimes.size(), binTimes.data(), binMeans.data(),
                                       xErr.data(), binErrors.data());
    g->SetLineWidth(2);
    g->SetMarkerStyle(21);
    g->GetXaxis()->SetTimeDisplay(1);
    g->SetTitle(Form("ARF; Time;Counts per hour (%dh average )", hoursBin));
    g->GetXaxis()->SetTimeFormat("%d/%m");
    g->GetXaxis()->SetTimeOffset(0,"gmt");
    return g;
}


void ARF() {
    struct DatasetConfig {
        std::string path;
        std::vector<std::pair<int,int>> rois;
        std::string comment;
    };

    std::vector<DatasetConfig> datasets = {
        {"/home/antoine/Documents/These/1-Radon/Radon_detector/3-Follow-Rn/ARF/250716_174300",{ {510,580},{590, 650},{770,840} },"ARF first use in July 2025"},
        {"/home/antoine/Documents/These/1-Radon/Radon_detector/3-Follow-Rn/ARF/251121_104230",{ {520,590},{600,670},{780,850} },"ARF verification issue, not enogh cold"},
        

    };

    std::vector<double> all_timestamp, roi1Rates, roi2Rates, roi3Rates;
    std::vector<double> roi1RateErr, roi2RateErr, roi3RateErr;
    std::vector<TH1F*> allHists;
    std::vector<TCanvas*> allCanvas;

    for(const auto& dataset : datasets) {
        std::string folderPath = dataset.path;
        std::string folderName = fs::path(folderPath).filename().string();
        std::cout<<"\n=== Processing dataset: "<<folderName<<" ===\n";

        std::vector<std::string> bkupFiles;
        for(const auto& entry : fs::directory_iterator(folderPath))
            if(entry.is_regular_file() && isValidBackupFile(entry.path().filename().string()))
                bkupFiles.push_back(entry.path().string());

        if(bkupFiles.empty()) { std::cout<<"No backup files found.\n"; continue; }

        std::sort(bkupFiles.begin(), bkupFiles.end(),
                  [](const std::string&a,const std::string&b){return extractBackupNumber(a)<extractBackupNumber(b);});

        time_t baseTimestamp = convertToTimestamp(folderName);
        if(baseTimestamp<0){ std::cerr<<"Invalid timestamp!\n"; continue; }

        std::vector<double> prevSpec(1024,0);
        double cumulativeTime=0, prevAcq=0;

        std::string totalFile="";
        for(const auto& entry: fs::directory_iterator(folderPath)) {
            if(entry.is_regular_file()) {
                std::string name = entry.path().filename().string();
                if(name.find("bkup")==std::string::npos && name.find(".iri")!=std::string::npos) {
                    totalFile = entry.path().string(); break;
                }
            }
        }

        if(!totalFile.empty()) {
            std::vector<double> totalSpec = readSpectrum(totalFile);
            TH1F* h = new TH1F(Form("h_%s",folderName.c_str()),
                               Form("Total Spectrum %s;Bin;Counts",folderName.c_str()),1024,0,1024);
            for(int b=0;b<1024;++b) h->SetBinContent(b+1,totalSpec[b]);
            allHists.push_back(h);

            TCanvas* c = new TCanvas(Form("c_%s",folderName.c_str()),
                                     Form("Spectrum %s",folderName.c_str()),800,600);
            allCanvas.push_back(c);
            c->cd();
            h->SetLineColor(kBlue); h->Draw("HIST"); c->Update();
        }

        for(size_t i=0;i<bkupFiles.size();++i) {
            const std::string& file = bkupFiles[i];
            double acq = readAcquisitionTime(file);
            double deltaT = (i==0)?acq:(acq-prevAcq);
            if(deltaT<=0){ std::cerr<<"Delta time error\n"; continue; }
            prevAcq = acq;

            auto curSpec = readSpectrum(file);
            std::vector<double> incrSpec(1024);
            if(i==0) incrSpec = curSpec;
            else for(int j=0;j<1024;++j) incrSpec[j] = curSpec[j]-prevSpec[j];

            auto roiRes = integrateROIs(incrSpec,dataset.rois,deltaT);
            time_t ts = baseTimestamp + static_cast<time_t>(cumulativeTime);
            all_timestamp.push_back(static_cast<double>(ts));
            cumulativeTime += deltaT;

            roi1Rates.push_back(roiRes["ROI_1"].rate); roi1RateErr.push_back(roiRes["ROI_1"].rateError);
            roi2Rates.push_back(roiRes["ROI_2"].rate); roi2RateErr.push_back(roiRes["ROI_2"].rateError);
            roi3Rates.push_back(roiRes["ROI_3"].rate); roi3RateErr.push_back(roiRes["ROI_3"].rateError);

            prevSpec = curSpec;
        }

        std::cout<<"DONE with folder: "<<folderName<<"\n";
    }

    if(!all_timestamp.empty()) {
        std::vector<double> xerr(all_timestamp.size(),3600/2);
        
        TCanvas* c1 = new TCanvas("c1","All ROI Rates",900,600);
        c1->SetLogy();
        c1->SetGridy();   
        gStyle->SetEndErrorSize(5);
        TGraphErrors* g2 = new TGraphErrors(all_timestamp.size(),all_timestamp.data(),roi2Rates.data(),
                                            xerr.data(),roi2RateErr.data());

        g2->SetLineColor(kMagenta); 
        g2->SetLineWidth(2); 
        g2->SetMarkerStyle(20);
        g2->SetMarkerSize(0);
        g2->Draw("AP");
        g2->GetYaxis()->SetRangeUser(1, 5000);

        TGraphErrors* g1 = new TGraphErrors(all_timestamp.size(),all_timestamp.data(),roi1Rates.data(),
                                            xerr.data(),roi1RateErr.data());
        g1->SetLineColor(kYellow+1); 
        g1->SetLineWidth(2); 
        g1->SetMarkerSize(0);
        g1->Draw("E same");

        TGraphErrors* g3 = new TGraphErrors(all_timestamp.size(),all_timestamp.data(),roi3Rates.data(),
                                            xerr.data(),roi3RateErr.data());
        g3->SetLineColor(kCyan+1); 
        g3->SetLineWidth(2); 
        g3->SetMarkerSize(0);
        g3->Draw("E same");

        g2->SetTitle("ARF; Time;Counts per hour");
        g2->GetXaxis()->SetTimeDisplay(1);
        g2->GetXaxis()->SetTimeFormat("%m/%d");

        TLatex latex;
        latex.SetNDC();          // Use normalized coordinates (0–1)
        latex.SetTextSize(0.04); // Adjust as needed

        latex.SetTextColor(kYellow+1);
        latex.DrawLatex(0.80, 0.85, "^{210}Po");

        latex.SetTextColor(kMagenta);
        latex.DrawLatex(0.80, 0.80, "^{218}Po");

        latex.SetTextColor(kCyan+1);
        latex.DrawLatex(0.80, 0.75, "^{214}Po");

        c1->Update();




        // --- Binned graph with errors ---
        long long tMin = static_cast<long long>(*std::min_element(all_timestamp.begin(),all_timestamp.end()));
        long long tMax = static_cast<long long>(*std::max_element(all_timestamp.begin(),all_timestamp.end()));
        int hoursBin = 3;
        

        TCanvas* c2 = new TCanvas("c2","Binned ROI Rates",900,600);
        c2->SetLogy();
        c2->SetGridy();   
        TGraphErrors* g11 = makeBinnedGraph(all_timestamp, roi1Rates, hoursBin, tMin, tMax);
        TGraphErrors* g22 = makeBinnedGraph(all_timestamp, roi2Rates, hoursBin, tMin, tMax);
        TGraphErrors* g33 = makeBinnedGraph(all_timestamp, roi3Rates, hoursBin, tMin, tMax);

        g11->SetLineColor(kYellow+1); 
        g11->SetMarkerColor(kYellow+1);
        g11->SetMarkerSize(0);

        g22->SetLineColor(kMagenta); 
        g22->SetMarkerColor(kMagenta);
        g22->SetMarkerSize(0);
        
        g33->SetLineColor(kCyan+1);
        g33->SetMarkerColor(kCyan+1);
        g33->SetMarkerSize(0);

        g22->Draw("AP");
        g22->GetXaxis()->SetTimeDisplay(1);
        g22->GetXaxis()->SetTimeFormat("%d/%m");
        g22->GetYaxis()->SetRangeUser(1, 5000);
        
        g11->Draw("P same");
        g33->Draw("P same");
        g11->SetTitle("Binned ROI Rates;Time;Counts per hour");

        
        latex.SetNDC();          // Use normalized coordinates (0–1)
        latex.SetTextSize(0.04); // Adjust as needed

        latex.SetTextColor(kYellow+1);
        latex.DrawLatex(0.80, 0.85, "^{210}Po");

        latex.SetTextColor(kMagenta);
        latex.DrawLatex(0.80, 0.80, "^{218}Po");

        latex.SetTextColor(kCyan+1);
        latex.DrawLatex(0.80, 0.75, "^{214}Po");
        c2->Update();


        // TLine ! 
        double f1, f2, f3, f4;
        
        


    } else std::cout<<"No points to plot\n";

    std::cout<<"\n=== ALL DATASETS FINISHED ===\n";
}
