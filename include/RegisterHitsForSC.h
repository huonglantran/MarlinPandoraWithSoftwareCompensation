/**
 *  @file   MarlinPandora/include/RegisterHitsForSC.h
 * 
 *  @brief  Header file for the register hits for SC algorithm class.
 * 
 *  $Log: $
 */
#ifndef REGISTER_HITS_FOR_SC_H
#define REGISTER_HITS_FOR_SC_H 1

#include "Pandora/Algorithm.h"

/**
 *  @brief  RegisterHitsForSC class
 */
class RegisterHitsForSC : public pandora::Algorithm
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
    RegisterHitsForSC();

    /**
      *  @breif Destructor
      */  
    ~RegisterHitsForSC();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ClusterType(const pandora::CaloHitList &caloHitList, int &isECalCluster, int &isHCalCluster) const;
    pandora::StatusCode ExtractCaloHits(const pandora::CaloHitList &caloHitList, std::vector<float> &ECalHitEnergies, std::vector<float> &HCalHitEnergies) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    std::string m_myRootFileName;
//    float m_HotHadCellEnergy;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *RegisterHitsForSC::Factory::CreateAlgorithm() const
{
    return new RegisterHitsForSC();
}

#endif //#ifndef REGISTER_HITS_FOR_SOFTWARECOMPENSATION_H
