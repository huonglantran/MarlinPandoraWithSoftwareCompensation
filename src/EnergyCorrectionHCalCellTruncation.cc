/**
 *  @file   MarlinPandoraSC/src/EnergyCorrectionHCalCellTruncation.cc
 * 
 *  @brief  Implementation of the HCalCellTruncation energy correction function algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "EnergyCorrectionHCalCellTruncation.h"
#include "PandoraPFANewProcessor.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

using namespace pandora;

EnergyCorrectionHCalCellTruncation::EnergyCorrectionHCalCellTruncation() :
    m_HCalCellTruncationEnergy(NULL)
{
  //fCluster = new TFile("ClusterInfo.root","recreate");
  //hClusterE = new TH1F("clusterEnergy", "clusterEnergy", 100, 0, 100);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionHCalCellTruncation::MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    if(NULL == pCluster)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if(0 == pCluster->GetNCaloHits())
    {
        correctedHadronicEnergy = 0.f;
        return STATUS_CODE_SUCCESS;
    }

    bool isHCalCluster(false), isECalCluster(false);
    
    //std::cout << "The hadronic energy for this cluster is:" << pCluster->GetHadronicEnergy() << std::endl;

    pandora::CaloHitList clusterCaloHitList;
    pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterCaloHitList);

    //float clusterEMenergy = pCluster->GetElectromagneticEnergy();
    float clusterHADenergy = pCluster->GetHadronicEnergy();
    //float clusterEnergy = clusterEMenergy + clusterHADenergy;

    //std::cout << "The EM energy for this cluster is:" << clusterEMenergy
    //<< " HAD energy " << clusterHADenergy
    //	      << std::endl;


    //std::cout << "clusterCaloHitList size " << clusterCaloHitList.size() << std::endl;

    this->clusterType(clusterCaloHitList, isECalCluster, isHCalCluster);

    if (isECalCluster)
        EnergyCorrectionHCalCellTruncation::ECalClusterEnergyCorrectionFunction(clusterCaloHitList, correctedHadronicEnergy);

    else if (isHCalCluster){
      EnergyCorrectionHCalCellTruncation::HCalCellTruncationClusterEnergyCorrectionFunction(clusterHADenergy,clusterCaloHitList, correctedHadronicEnergy);
      //hClusterE->Fill(clusterHADenergy);
    }
    else
      EnergyCorrectionHCalCellTruncation::HCalCellTruncationHCalECalSplitClusterEnergyCorrectionFunction(clusterHADenergy,clusterCaloHitList, correctedHadronicEnergy);
    
    return STATUS_CODE_SUCCESS;

}


//------------------------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionHCalCellTruncation::clusterType(const pandora::CaloHitList &caloHitList, bool &isECalCluster, bool &isHCalCluster) const
{
    int nECalHits(0), nHCalHits(0);

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;

        if (HCAL == pCaloHit->GetHitType())
            nHCalHits++;

        else if (ECAL == pCaloHit->GetHitType())
            nECalHits++;

        //std::cout << "pCaloHit->GetHitType():   " << pCaloHit->GetHitType() << std::endl;
        //std::cout << "nECalHits:                " << nECalHits << std::endl;
        //std::cout << "nHCalHits:                " << nHCalHits << std::endl;
    }

    if (nECalHits != 0 && nHCalHits == 0)
        isECalCluster = true;

    else if (nHCalHits != 0 && nECalHits == 0)
        isHCalCluster = true;

    //std::cout << "nECalHits:    " << nECalHits << std::endl;
    //std::cout << "nHCalHits:    " << nHCalHits << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionHCalCellTruncation::ECalClusterEnergyCorrectionFunction(const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
    //std::cout << "========ECAL CONTAINED CLUSTER========" << std::endl;

    for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        energyCorrection += pCaloHit->GetHadronicEnergy();
        //std::cout << "Calo Hit Energy: " << pCaloHit->GetInputEnergy() << std::endl;
    }

    //std::cout << "energyCorrection:   " << energyCorrection <<std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionHCalCellTruncation::HCalCellTruncationClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{
  float E_HCalCellTruncation(0.f);

  for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
        const pandora::CaloHit *pCaloHit = *iter;
        //std::cout << "Calo Hit Energy: " << pCaloHit->GetInputEnergy() << std::endl;
        //std::cout << "Calo Hit Type:   " << pCaloHit->GetHitType() << std::endl;

	float hitEnergy = pCaloHit->GetHadronicEnergy();
	//float hitEnergy = pCaloHit->GetInputEnergy();
	
	if (hitEnergy>m_HCalCellTruncationEnergy.at(0)) 
	  hitEnergy = m_HCalCellTruncationEnergy.at(0);

	E_HCalCellTruncation += hitEnergy;
    }

    energyCorrection = E_HCalCellTruncation;

    //std::cout << "energyCorrection:     " << energyCorrection <<std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionHCalCellTruncation::HCalCellTruncationHCalECalSplitClusterEnergyCorrectionFunction(float clusterEnergyEstimation, const pandora::CaloHitList &caloHitList, float &energyCorrection) const
{

  float E_HCalCellTruncation(0.f);
  float E_EM(0.f);
  float E_HAD(0.f);
  float E_HAD_HCalCellTruncation(0.f);

  //std::cout << "Energy correction for split clusters: clusterEnergyEstimation " << clusterEnergyEstimation << std::endl;

  for(pandora::CaloHitList::const_iterator iter = caloHitList.begin() , endIter = caloHitList.end() ; endIter != iter ; ++iter)
    {
      const pandora::CaloHit *pCaloHit = *iter;
      //std::cout << "Calo Hit Energy: " << pCaloHit->GetInputEnergy() << std::endl;
      //std::cout << "Calo Hit Type:   " << pCaloHit->GetHitType() << std::endl;

	//If HCAL: correct
        if (HCAL == pCaloHit->GetHitType())
        {
	  float hitEnergy = pCaloHit->GetHadronicEnergy();
	  //float hitEnergy_Had = pCaloHit->GetHadronicEnergy();
	  
	  if (hitEnergy>m_HCalCellTruncationEnergy.at(0)) 
	    hitEnergy = m_HCalCellTruncationEnergy.at(0);
	  
	  E_HAD += hitEnergy;
	  E_HAD_HCalCellTruncation += hitEnergy;
	  E_HCalCellTruncation += hitEnergy;
        }	
        else
	  {           
	  E_EM += pCaloHit->GetHadronicEnergy();
	  E_HCalCellTruncation += pCaloHit->GetHadronicEnergy();
        }
    }

    energyCorrection = E_HCalCellTruncation;

    //std::cout << "The corrected hadronic energy for a cluster of threshold hits:" << std::endl;
    //std::cout << "E_EM = " << E_EM << " E_HAD = " << E_HAD << " E_HAD_HCalCellTruncation " << E_HAD_HCalCellTruncation << std::endl;
    //std::cout << "energyCorrection:     " << energyCorrection << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyCorrectionHCalCellTruncation::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "HCalCellTruncationEnergy", m_HCalCellTruncationEnergy));

    //std::cout << "====== Read Settings =====" << std::endl;
    std::cout << "m_HCalCellTruncationEnergy : " << m_HCalCellTruncationEnergy.at(0) << std::endl;

    return STATUS_CODE_SUCCESS;
}

//I =1; J = 2, K = 12
//I = 1; j = 15; k = 12

