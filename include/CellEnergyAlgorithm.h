/**
 *  @file   MarlinPandora/include/CellEnergyAlgorithm.h
 * 
 *  @brief  Header file for the cell energy algorithm class.
 * 
 *  $Log: $
 */
#ifndef CELL_ENERGY_ALGORITHM_H
#define CELL_ENERGY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

//#include "TFile.h"
//#include "TH1F.h"
//#include "TTree.h"

/**
 *  @brief  CellEnergyAlgorithm class
 */
class CellEnergyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
      *  @brief  Default constructor
      */
    CellEnergyAlgorithm();

    /**
      *  @breif Destructor
      */  
    ~CellEnergyAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ClusterType(const pandora::CaloHitList &caloHitList, int &isECalCluster, int &isHCalCluster) const;
    pandora::StatusCode Display(const pandora::CaloHitList &caloHitList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode ClusterDistributionMetric(const pandora::CaloHitList &caloHitList, float &clusterDistributionMetric, float &clusterHadronicEnergy, int &nHitsAbvMean) const;

    // Member variables here
    std::string m_RootFileName;
//    float m_HotHadCellEnergy;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CellEnergyAlgorithm::Factory::CreateAlgorithm() const
{
    return new CellEnergyAlgorithm();
}

#endif // #ifndef CELL_ENERGY_ALGORITHM_H
