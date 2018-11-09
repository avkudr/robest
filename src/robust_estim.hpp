/**
 *  @brief Robust estimation library
 *
 *  @author  Andrey Kudryavtsev (avkudr.github.io)
 *  @date    01/03/2018
 *  @version 1.0
 */

#ifndef ROBUST_ESTIMATOR_H
#define ROBUST_ESTIMATOR_H

#include <vector>
#include <random>
#include <numeric>
#include <functional>
#include <algorithm>
#include <iostream>
#include <iterator>

namespace robest
{

class EstimationProblem
{

  public:
    // Functions to overload in your class:
    virtual double estimErrorForSample(int i) = 0;
    virtual void estimModelFromSamples(std::vector<int> samplesIdx) = 0;
    virtual int getTotalNbSamples() const = 0;

    int getNbParams() const { return nbParams; }
    int getNbMinSamples() const { return nbMinSamples; }

  protected:
    void setNbParams(int i) { nbParams = i; }
    void setNbMinSamples(int i) { nbMinSamples = i; }

  private:
    int nbParams = -1;
    int nbMinSamples = -1;
};

static std::random_device rd;  // random device engine, usually based on /dev/random on UNIX-like systems
static std::mt19937 rng(rd()); // initialize Mersennes' twister using rd to generate the seed
//static std::mt19937 rng((unsigned int) - 1);

// Base class for Robust Estimators
class AbstractEstimator
{
  public:
    // Generate X !different! random numbers from 0 to N-1
    // X - minimal number of samples
    // N - total number of samples
    // generated numbers are indices of data points

    std::vector<int> randomSampleIdx()
    {
        int minNbSamples = problem->getNbMinSamples();
        int totalNbSamples = problem->getTotalNbSamples();

        std::vector<int> allIdx(totalNbSamples);
        std::iota(std::begin(allIdx), std::end(allIdx), 0); // Fill with 0, 1, ..., totalNbSamples-1.

        // shuffle the elements of allIdx : Fisher–Yates shuffle
        for (int i = 0; i < minNbSamples; i++)
        {
            std::uniform_int_distribution<int> dist(0, totalNbSamples - i - 1);
            int randInt = dist(rng);
            std::swap(allIdx[totalNbSamples - i - 1], allIdx[randInt]);
        }

        //take last <minNbSamples> elements
        std::vector<int> idx(allIdx.end() - minNbSamples, allIdx.end());
        //        std::cout << "[";
        //        for (auto i : idx) std::cout << i << ", ";
        //        std::cout << "]\n";

        return idx;
    }

    int estimOfIter(int dataSize, float alpha = 0.99, float gamma = 0.60)
    {
        return (std::log(alpha)) / (std::log(1 - std::pow(gamma, dataSize)));
    }

    double getInliersFraction() const { return inliersFraction; }
    std::vector<int> getInliersIndices() const { return inliersIdx; }

  protected:
    void getInliers(double thres)
    {
        int totalNbSamples = problem->getTotalNbSamples();
        inliersIdx.clear();
        for (int j = 0; j < totalNbSamples; j++)
        {
            double error = problem->estimErrorForSample(j);
            
            if (error * error < thres * thres)
            {
                inliersIdx.push_back(j);
            }
        }
        this->inliersFraction = (double)(inliersIdx.size()) / (double)(totalNbSamples);
    }

    EstimationProblem *problem;
    std::vector<int> bestIdxSet;
    std::vector<int> inliersIdx;
    double inliersFraction = -1.0;

};

class RANSAC : public AbstractEstimator
{
  public:
    RANSAC()
    {
    }

    void solve(EstimationProblem *pb, double thres = 0.1, int nbIter = 10000)
    {
        problem = pb;

        int totalNbSamples = problem->getTotalNbSamples();
        for (int i = 0; i < nbIter; i++)
        {
            std::vector<int> indices = randomSampleIdx();

            problem->estimModelFromSamples(indices);

            //getInliersNb
            int nbInliers = 0;
            for (int j = 0; j < totalNbSamples; j++)
            {
                double error = problem->estimErrorForSample(j);
                error = error * error;
                if (std::fabs(error) < thres)
                {
                    nbInliers++;
                }
            }

            double inliersFraction = (double)(nbInliers) / (double)(problem->getTotalNbSamples());
            if (inliersFraction > this->inliersFraction)
            {
                this->inliersFraction = inliersFraction;
                this->bestIdxSet = indices;
            }
        }
        problem->estimModelFromSamples(bestIdxSet);
        getInliers(thres);
        problem->estimModelFromSamples(inliersIdx);
    }
};

class MSAC : public AbstractEstimator
{
  public:
    MSAC()
    {
        sumSqErr = 1000000000.0;
    }

    void solve(EstimationProblem *pb, double thres = 0.1, int nbIter = 10000)
    {
        problem = pb;

        int totalNbSamples = problem->getTotalNbSamples();

        if (nbIter == 0)
        {
            nbIter = this->estimOfIter(totalNbSamples);
            std::cout << std::endl;
            std::cout << nbIter << std::endl;
        }

        for (int i = 0; i < nbIter; i++)
        {
            std::vector<int> indices = randomSampleIdx();
            problem->estimModelFromSamples(indices);

            //getInliersNb
            double sumSqErr = 0;
            int nbInliers = 0;
            for (int j = 0; j < totalNbSamples; j++)
            {
                double error = problem->estimErrorForSample(j);
                if (error * error < thres * thres)
                {
                    sumSqErr += error * error;
                    nbInliers++;
                }
                else
                {
                    sumSqErr += thres * thres;
                }
            }

            double inliersFraction = (double)(nbInliers) / (double)(problem->getTotalNbSamples());
            //if (sumSqErr < this->sumSqErr && inliersFraction > this->inliersFraction)
            if (sumSqErr < this->sumSqErr )
            {
                this->inliersFraction = inliersFraction;
                this->sumSqErr = sumSqErr;
                this->bestIdxSet = indices;
            }
        }
        
        problem->estimModelFromSamples(bestIdxSet);
        getInliers(thres);
        //problem->estimModelFromSamples(inliersIdx);
    }
  private:
    double sumSqErr;
};
/*
class MLESAC : public AbstractEstimator
{
  public:
    MLESAC()
    {
        
    }

    void solve(EstimationProblem *pb, double thres = 0.1, int nbIter = 10000, double Sigma, int nEMIter)
    {
        problem = pb;

        int totalNbSamples = problem->getTotalNbSamples();

        if (nbIter == 0)
        {
            nbIter = this->estimOfIter(totalNbSamples);
            std::cout << std::endl;
            std::cout << nbIter << std::endl;
        }

        for (int i = 0; i < nbIter; i++)
        {
            std::vector<int> indices = randomSampleIdx();
            problem->estimModelFromSamples(indices);

            for (int j = 0; j < totalNbSamples; j++)
            {
                double error = problem->estimErrorForSample(j);
                
                dMix = 0.5;
                
                for (int  j = 0; j > nEMIter; j++)
                {
                    dErrorInlierProb = dMix * exp(-std::pow(error,2)/(2*std::pow(Sigma,2)))/(Sigma * std::sqrt(2 * 3.1415));
                    dErrorOutlierProb = (1 - dMix) / v;
                }
            }

            double inliersFraction = (double)(nbInliers) / (double)(problem->getTotalNbSamples());
            //if (sumSqErr < this->sumSqErr && inliersFraction > this->inliersFraction)
            if (sumSqErr < this->sumSqErr )
            {
                this->inliersFraction = inliersFraction;
                this->sumSqErr = sumSqErr;
                this->bestIdxSet = indices;
            }
        }
        
        problem->estimModelFromSamples(bestIdxSet);
        getInliers(thres);
        //problem->estimModelFromSamples(inliersIdx);
    }
  private:
    double dMix;
    double dErrorInlierProb;
    double dErrorOutlierProb;
    double dInlierProb;
    double v;
};
*/
class LMedS : public AbstractEstimator
{
  public:
    LMedS()
    {
        med = 1000000.0;
    }

    template <typename T>
    double median(std::vector<T> &v)
    {
        std::sort(v.begin(), v.end());
        if (v.size() % 2 == 0)
        {
            return (double)(v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
        }
        else
        {
            return v[v.size() / 2];
        }
    }

    void solve(EstimationProblem *pb, double thres = 0.1, int nbIter = 1000)
    {
        problem = pb;
        int totalNbSamples = problem->getTotalNbSamples();

        for (int i = 0; i < nbIter; i++)
        {
            std::vector<int> indices = randomSampleIdx();

            problem->estimModelFromSamples(indices);

            std::vector<double> errorsVec(totalNbSamples);
            for (int j = 0; j < problem->getTotalNbSamples(); j++)
            {
                double error = problem->estimErrorForSample(j);
                //errorsVec[j] = error; // error must be squared!
                errorsVec[j] = error * error; // error must be squared!
            }
            double med = median(errorsVec);
            if (med < this->med)
            {
                this->med = med;
                this->bestIdxSet = indices;
            }
        }

        problem->estimModelFromSamples(bestIdxSet);
        getInliers(thres);
        problem->estimModelFromSamples(inliersIdx);
    }

  private:
    double med;
};

} // namespace robest
#endif // ROBUST_ESTIMATOR_H
