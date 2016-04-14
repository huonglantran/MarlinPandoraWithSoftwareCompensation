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

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Is the cluster contained within the ECal or HCal or split?
     */
    pandora::StatusCode clusterType(const pandora::CaloHitList &caloHitList, bool &isECalCluster, bool &isHCalCluster) const;

    /**
     *  @brief  To find hit energy density
     */
    pandora::StatusCode FindDensity(const pandora::CaloHit *pCaloHit, float &hitEnergyDensity) const;

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

    pandora::FloatVector    m_SCEnergyConstants1;                 //
    pandora::FloatVector    m_SCEnergyConstants2;                 //
    bool                    m_cutClusterHadEnergy;                // Flag to apply cluster energy cut
    float                   m_clusterMinHadEnergy;                // Minimum hadronic energy for a cluster to apply software compensation to  
    bool                    m_cutNumberOfHitsInCluster;           // Flag to apply number of hits cut
    int                     m_minNumberOfHitsInCluster;           // Minimum number of hits for a cluster to apply software compensation
    bool                    m_applyOnlyToTrackAssociatedClusters; // Flag to turn on when desire to apply Software Compensation only to cluster with at least one track associated
    bool                    m_cutCglobal;                         // Flag to apply cut on cluster energy distribution variable
    float                   m_maxCglobal;                         // Maximum value of cluster energy distribution variable to apply software compensation
    bool                    m_cutClusterTopologyVariable;         // Flag to apply cluster topology variable cut
    float                   m_minClusterTopologyVariable;         // Minimum value of cluster topology variable to apply software compensation

    bool                    m_cheating;                           // Flag for cheating 
    float                   m_trueEnergy;                         // Cheating: for single particle samples one can give true particle energy

};

#endif // #ifndef SC_ENERGY_CORRECTION_PLUGINS_H
