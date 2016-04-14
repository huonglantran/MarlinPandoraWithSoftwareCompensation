/**
 *  @file   MarlinPandoraSC/include/EnergyCorrectionSC.h
 * 
 *  @brief  Header file for the SC energy correction plugin algorithm class.
 * 
 *  $Log: $
 */

#ifndef HCALCELLTRUNCATION_ENERGY_CORRECTION_PLUGINS_H 
#define HCALCELLTRUNCATION_ENERGY_CORRECTION_PLUGINS_H 1

#include "Plugins/EnergyCorrectionsPlugin.h"

class TFile;
class TH1F;
class TTree;

/**
 *   @brief  EnergyCorrectionHCalCellTruncation class. 
 */

class EnergyCorrectionHCalCellTruncation : public pandora::EnergyCorrectionPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    EnergyCorrectionHCalCellTruncation();

    /**
     *  @brief  Algorithm to implement HCalCellTruncation hadronic energy calculation
     */
    pandora::StatusCode MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const;

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
    pandora::StatusCode HCalCellTruncationClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const;

    /**
     *  @brief  Calculation of the hadronic energy, if cluster split between in ECal and HCal, 
     *          based on number of calo hits in HCal and raw hadronic energy in the ECal
     */
    pandora::StatusCode HCalCellTruncationHCalECalSplitClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const;

typedef std::vector<float> FloatVec;

    FloatVec         m_HCalCellTruncationEnergy;            //

    //Lan add to make cluster energy distribution
    //TFile *fCluster;
    //TH1F *hClusterE;

};

#endif // #ifndef HCALCELLTRUNCATION_ENERGY_CORRECTION_PLUGINS_H
