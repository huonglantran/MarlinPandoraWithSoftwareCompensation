/**
 *  @file   MarlinPandoraSC/include/EnergyCorrectionSC.h
 * 
 *  @brief  Header file for the SC energy correction plugin algorithm class.
 * 
 *  $Log: $
 */

#ifndef SC_ENERGY_CORRECTION_PLUGINS_H 
#define SC_ENERGY_CORRECTION_PLUGINS_H 1

#include "Plugins/EnergyCorrectionsPlugin.h"

class TFile;
class TH1F;
class TTree;

/**
 *   @brief  EnergyCorrectionSC class. 
 */

class EnergyCorrectionSC : public pandora::EnergyCorrectionPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    EnergyCorrectionSC();

    /**
     *  @brief  Algorithm to implement SC hadronic energy calculation
     */
    pandora::StatusCode MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const;

    //float FindDensity(float hitEnergy);
    
private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Is the cluster contained within the ECal or HCal or split?
     */
    pandora::StatusCode clusterType(const pandora::CaloHitList &caloHitList, bool &isECalCluster, bool &isHCalCluster) const;

    /**
     *  @brief  Calculation of the hadronic energy, if cluster contained in ECal.  Returned value is sum of 
     *          raw hadronic energy in the ECal
     */
    pandora::StatusCode ECalClusterEnergyCorrectionFunction(const pandora::CaloHitList &caloHitList, float &energyCorrection) const;

    /**
     *  @brief  Calculation of the hadronic energy, if cluster contained in HCal, based on number of calo hits
     */
    pandora::StatusCode SCClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const;

    /**
     *  @brief  Calculation of the hadronic energy, if cluster split between in ECal and HCal, 
     *          based on number of calo hits in HCal and raw hadronic energy in the ECal
     */
    pandora::StatusCode SCHCalECalSplitClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const;

    /**
     *  @brief  Calculate energy density in units of MIP per cell
     */
    pandora::StatusCode FindDensity(const pandora::CaloHit *const pCaloHit, float &energyDensity) const;

    /**
     *  @brief  ... 
     */
    pandora::StatusCode CleanCluster(const pandora::Cluster *const pCluster, float &correctedHadronicEnergy) const;
    float GetHadronicEnergyInLayer(const pandora::OrderedCaloHitList &orderedCaloHitList, const unsigned int pseudoLayer) const;

    pandora::FloatVector          m_SCEnergyConstants1;            //
    pandora::FloatVector          m_SCEnergyConstants2;            //
    bool                          m_cheating;  
    float                         m_trueEnergy;
    float                         m_clusterMinHadEnergy;           // Minimum hadronic energy for a cluster to apply software compensation to
    float           m_minCleanHitEnergy;                ///< Min calo hit hadronic energy to consider cleaning hit/cluster
    float           m_minCleanHitEnergyFraction;        ///< Min fraction of cluster energy represented by hit to consider cleaning
    float           m_minCleanCorrectedHitEnergy;       ///< Min value of new hit hadronic energy estimate after cleaning
};

#endif // #ifndef SC_ENERGY_CORRECTION_PLUGINS_H
